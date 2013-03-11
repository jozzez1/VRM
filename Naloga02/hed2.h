// hed2.h
////////////

#ifndef __HEADER_VRM2
#define __HEADER_VRM2

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_gamma.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <stdlib.h>

typedef struct
{
	double h,     // space discretization
	       L,     // lambda ... parameter of potential
	       a,     // f(x) -> f(x - a)
	       t;     // time step

	int N,        // rank of the matrices
	    tryrac,   // accuracy "boolean"
	    p,        // approximation order of the 2nd derivative
	    n,        // which vector do we want for initial condition
	    d,        // animation flag
	    k,        // current time iteration
	    T;        // number of time iterations

	double * H,   // Hamiltonian, in the end matrix v
	       * v,   // eigenvector matrix in Lanczos space
	       * w,   // Lanczos transformation matrix
	       * DT,  // dominant of tridiagonal
	       * ST,  // subdiagonal of tridiagonal
	       * E;   // eigenvalue array -- energies ...
	
	char * dat;
} hod;

// potential
double pot (hod * u, int i)
{
	double x = u->h * (i - 0.5*u->N),
	       V = 0.5*x*x + u->L*pow (x, 4);

	return V;
}

// tridiagonal hamiltonian, O(h^2) precision
void hamilton2 (hod * u)
{
	int i;

	gsl_matrix_view m = gsl_matrix_view_array (u->H, u->N, u->N);

	double a = 1.0/pow(u->h, 2),
	       b = (-0.5)*a;

	for (i = 0; i <= u->N-2; i++)
		gsl_matrix_set (&m.matrix, i, i+1, b);

	for (i = 1; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-1, b);

	for (i = 0; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i, a + pot(u, i));

}

// five-diagonal hamiltonian, O(h^4) precision
void hamilton4 (hod * u)
{
	int i;

	gsl_matrix_view m = gsl_matrix_view_array (u->H, u->N, u->N);

	double a = 1.0/pow(u->h, 2),
	       b2 = (1.0/24)*a,
	       b1 = (-2.0/3)*a,
	       d = 1.25*a;

	for (i = 0; i <= u->N-3; i++)
		gsl_matrix_set (&m.matrix, i, i+2, b2);

	for (i = 0; i <= u->N-2; i++)
		gsl_matrix_set (&m.matrix, i, i+1, b1);

	for (i = 0; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i, d + pot(u, i));

	for (i = 1; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-1, b1);

	for (i = 2; i <= u->N-2; i++)
		gsl_matrix_set (&m.matrix, i, i-2, b2);
}

// seven-diagonal approximation, O(h^6) precision
void hamilton6 (hod * u)
{
	int i;

	gsl_matrix_view m = gsl_matrix_view_array (u->H, u->N, u->N);

	double a = 1.0/pow (u->h, 2),
	       b3 = (-1)*(0.5/90)*a,
	       b2 = (0.5)*(3.0/20)*a,
	       b1 = (-0.5)*1.5*a,
	       d = 0.5*(49.0/18)*a;

	for (i = 0; i <= u->N-4; i++)
		gsl_matrix_set (&m.matrix, i, i+3, b3);

	for (i = 0; i <= u->N-3; i++)
		gsl_matrix_set (&m.matrix, i, i+2, b2);

	for (i = 0; i <= u->N-2; i++)
		gsl_matrix_set (&m.matrix, i, i+1, b1);

	for (i = 0; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i, d + pot(u, i));

	for (i = 1; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-1, b1);

	for (i = 2; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-2, b2);

	for (i = 3; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-3, b3);
}

// hamilton, nine-diagonal approximation, O(h^8) precision
void hamilton8 (hod * u)
{
	int i;

	gsl_matrix_view m = gsl_matrix_view_array (u->H, u->N, u->N);

	double a = 1.0/pow (u->h, 2),
	       b4 = (0.5/560)*a,
	       b3 = (-0.5)*(8.0/315)*a,
	       b2 = 0.1*a,
	       b1 = (-0.5)*1.6*a,
	       d = 0.5*(205.0/72)*a;

	for (i = 0; i <= u->N-5; i++)
		gsl_matrix_set (&m.matrix, i, i+4, b4);

	for (i = 0; i <= u->N-4; i++)
		gsl_matrix_set (&m.matrix, i, i+3, b3);

	for (i = 0; i <= u->N-3; i++)
		gsl_matrix_set (&m.matrix, i, i+2, b2);

	for (i = 0; i <= u->N-2; i++)
		gsl_matrix_set (&m.matrix, i, i+1, b1);

	for (i = 0; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i, d + pot(u, i));

	for (i = 1; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-1, b1);

	for (i = 2; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-2, b2);

	for (i = 3; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-3, b3);

	for (i = 4; i <= u->N-1; i++)
		gsl_matrix_set (&m.matrix, i, i-4, b4);
}

void swap (double * v1, double * v2, int N)
{
	int i;
	double x;
	for (i = 0; i <= N-1; i++)
	{
		x = v1[i];
		v1[i] = v2[i];
		v2[i] = x;
	}
}

double scalar (double * v1, double * v2, int N)
{
	int i;
	double re = 0;
	
	for (i = 0; i <= N-1; i++)
		re += v1[i] * v2[i];

	return re;
}

// Lanczos algorithm
void Lanczos (hod * u)
{
	int i;

	gsl_matrix_view m = gsl_matrix_view_array (u->w, u->N, u->N);

	double * b = (double *) malloc (u->N * sizeof (double)),
	       * a = (double *) malloc (u->N * sizeof (double)),
	       * v0 = (double *) calloc (u->N, sizeof (double)),
	       * v1 = (double *) calloc (u->N, sizeof (double)),
	       * w = (double *) calloc (u->N, sizeof (double));

	b[0] = 0;
	v1[0] = 1;

	for (i = 0; i <= u->N - 2; i++)
		gsl_matrix_set (&m.matrix, i, 0, v1[i]);

	for (i = 0; i <= u->N-2; i++)
	{
		cblas_dsymv (CblasRowMajor, CblasUpper, u->N,
				1.0, u->H, u->N, v1, 1, 0.0, w, 1);
		a[i] = scalar (w, v1, u->N);

		int j;
		for (j = 0; j <= u->N - 1; j++)
			w[j] -= (a[i]*v1[j] + b[i]*v0[j]);

		b[i+1] = sqrt (scalar (w, w, u->N));

		swap (v1, v0, u->N);

		for (j = 0; j <= u->N-1; j++)
		{
			v1[j] = w[j]/b[i+1];
			gsl_matrix_set (&m.matrix, j, i, v1[j]);
		}
	}

	cblas_dsymv (CblasRowMajor, CblasUpper, u->N,
			1.0, u->H, u->N, v1, 1, 0.0, w, 1);
	a[u->N-1] = scalar (w, v1, u->N);

	// we feed those parameters to the struct
	for (i = 0; i <= u->N - 2; i++)
	{
		u->DT [i] = a[i];
		u->ST [i] = b[i+1];
	}

	u->DT [u->N-1] = a[i];

	free (b);
	free (a);
	free (v0);
	free (v1);
	free (w);
}

// struct initializer
void init (hod * u, double h, double L, double a, double t, int N,
		int p, int d, int tryrac, int n, int T, char * dat)
{
	u->h = h;
	u->L = L;
	u->a = a;
	u->t = t;
	u->N = N;
	u->p = p;
	u->n = n;
	u->T = T;
	u->d = d;
	u->k = 0;
	u->tryrac = tryrac;

	int M = u->N * u->N,
	    D = u->N - 1;

	u->H = (double *) calloc (M, sizeof(double));
	u->v = (double *) malloc (M * sizeof (double));
	u->w = (double *) malloc (M * sizeof (double));
	u->DT = (double *) malloc (u->N * sizeof (double));
	u->ST = (double *) malloc (D * sizeof (double));
	u->E = (double *) malloc (u->N * sizeof (double));
	u->dat = (char *) malloc (15 * sizeof (char));

	if (dat == NULL)
		sprintf (u->dat, "nal-N%d-L%d-a%d",
				u->N, (int) (100 * u->L), (int) u->a);

	else
		strcpy (u->dat, dat);

	// we fill the hamiltonian
	switch (u->p)
	{
		case 2:
			hamilton2 (u);
			break;

		case 4:
			hamilton4 (u);
			break;

		case 6:
			hamilton6 (u);
			break;

		case 8:
			hamilton8 (u);
			break;

		default:
			printf ("Invalid option!\n");
			printf ("Selecting the five-diagonal option instead.\n");
			hamilton4 (u);
			break;
	}

	// we use the Lanczos method for hermitian matrices
	if (u->p != 2)
		Lanczos (u);
}

// Elenaki ise i jineki mou
void destroy (hod * u)
{
	free (u->H);
	free (u->v);
	free (u->w);
	free (u->DT);
	free (u->ST);
	free (u->E);

	free (u);
}


// now to diagonalize this baby with this handy f77 wrapper ...
// compile with: -llapack -lblas -lcblas -latlas
// for fortran code ALL arguments must be passed as pointers
static long MRRR (hod * u)
{
	// automatically tries for 8 decimal accuracy
	extern void dstemr_ (char * jobz, // calculation mode
			char * range,     // range of eigenvalues to be computed
			int * N,          // rank of the matrix to be computed
			double * D,       // (in/out) leading diagonal of the matrix
			double * E,       // (in/out) side diagonals of the matrix
			double * vl,      // lowest window for the eigenvalue
			double * vu,      // upper bound for the eigen value
			int * Il,         // again some bounds ... set to 0
			int * Iu,         // again bounds, depending on RANGE variable
			int * M,          // (out) number of eigenvalues found
			double * w,       // (out) vector of the found eigenvalues
			double * z,       // (out) array of eigenvectors -- (in the columns)
			int * ldz,        // leading dimension of the array
			int * nzc,        // number of eigenvectors to be held in z
			int * isuppz,     // (out) integer array ... dimension 2*max(1,M)
			int * tryrac,     // (in/out) logical (int or long int) for high accuracy
			double * work,    // (out) array of dimension lwork, which is the output
			int * lwork,      // (in) for JOBZ='V' has to be >= max (1, 12*N)
			double * iwork,   // (out) same as work
			int * liwork,     // (in) for JOBZ='V' has to be >= max (1, 10*N)
			int * Info);      // (out) =0 for successfull operation

	char JOBZ   =   'V',
	     RANGE  =   'A';

	int val = u->N-1,
	    i;

	double * d = (double *) malloc (u->N * sizeof(double)),
	       * e = (double *) malloc (val * sizeof (double));
	
	for (i = 0; i <= u->N - 2; i++)
	{
		d[i] = u->DT[i];
		e[i] = u->ST[i];
	}
	d[u->N - 1] = u->DT[i];

	double VL = 0,
	       VU = 1;

	int IL = 0,
	    IU = 1,
	    ispz = 2 * u->N,
	    M;

	// for 'w' I will use 'u->e' -- already an array
	// for 'z' I will use 'u->v' -- already an array
	
	int LDZ    = u->N,
	    NZC    = u->N,
	    TRYRAC = u->tryrac, // I hope this means 'yes, check for precision :D'
	    LWORK  = 18*u->N,
	    LIWORK = 10*u->N,
	    INFO;
	
	int * ISUPPZ = (int *) malloc (ispz * sizeof (int));

	double * WORK  = (double *) malloc (LWORK * sizeof (double)),
	       * IWORK = (double *) malloc (LIWORK * sizeof (double));

	dstemr_ (&JOBZ, &RANGE, &u->N, d, e, &VL, &VU, &IL, &IU, &M, u->E, u->v, &LDZ, &NZC,
			ISUPPZ, &TRYRAC, WORK, &LWORK, IWORK, &LIWORK, &INFO);

// if I deallocate but one of them I get a SEGSIGV ...
// until I figure what to do, I'll comment these out
//	free (d);
//	free (e);
//	free (ISUPPZ);
//	free (WORK);
//	free (IWORK);

	return INFO;
}

// calculate eigenvectors in the basis of Harmonic potential
void rotate (hod * u)
{
	// we have to multiply Lanczos and Eigenmatrix
	// H = w*v
	// w -- lanczos         -- it's ok ... row major ordering ;)
	// v -- diagonalization	-- because of fortran it is in column major ordering
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans, u->N, u->N, u->N,
			1.0, u->w, u->N, u->v, u->N, 0.0, u->H, u->N);

	// now the former hamiltonian matrix contains the eigenvectors
	// on its columns, and thus corresponds to the transposed matrix
	// of cannonical diagonalization.
}

void eigen_dump (hod *u)
{
		char * eigen = (char *) malloc (12 * sizeof (char));
		sprintf (eigen, "Energies-N%d.txt", u->N);

		FILE * fout = fopen (eigen, "w");
		free (eigen);

		int i;
		for (i = 0; i <= u->N - 1; i++)
			fprintf (fout, "% 4d % 15.8lf\n", i, u->E[i]);

		fclose (fout);
}

// diagonalize O(n^2)
void diag_MRRR (hod * u)
{
	int info = MRRR (u),
	    i;

	if (info != 0)
	{
		fprintf (stderr, "Error %d! Unable to complete.\n", info);
		exit (EXIT_FAILURE);
	}

	for (i = 0; i <= u->N - 1; i++)
		printf ("% 4d % 15.8lf\n", u->N-1-i, u->E[u->N-1-i]);

	printf ("Accuracy: %d\n", u->tryrac);

	// we just rotate to the correct base
	rotate (u);
}

// eigenvectors for harmonic potential
double phi (int n, double x)
{
	double norm = 1.0/(pow (M_PI, 0.25) * sqrt(gsl_sf_fact (n) * pow(2, n))),
	       S = 0;

	int i, N,
	    k = n%2;

	if (k == 0)
	{
		N = n/2;
		for (i = 0; i <= N; i++)
			S += pow((-1), N - i)*pow(2*x, 2*i)/(gsl_sf_fact(2*i)*gsl_sf_fact(N - i));

		S *= norm * gsl_sf_fact (n) * exp ((-0.5)*x*x);
	}

	else if (k == 1)
	{
		N = (n - 1)/2;
		for (i = 0; i <= N; i++)
			S += pow((-1), N - i)*pow(2*x, 2*i+1)/(gsl_sf_fact(2*i+1)*gsl_sf_fact(N-i));

		S *= norm * gsl_sf_fact (n) * exp ((-0.5)*x*x);
	}

	return S;
}

// we can put initial vector in unused u->DT and u->ST
// to propagate them through time
void init_v (hod * u)
{
	int i;

	// we initialize the real part
	for (i = 0; i <= u->N - 1; i++)
		u->DT [i] = phi (u->n, u->h * (i - 0.5*u->N) - u->a);

	// we have to modify the imaginary part -- ST
	u->ST = (double *) realloc ((double *) u->ST, u->N * sizeof (double));
	for (i = 0; i <= u->N - 1; i++)
		u->ST [i] = 0.0;

	// we shall also reuse the v and w:
	// v shall be DT in psi space
	// w shall be ST in psi space
	u->v = (double *) realloc ((double *) u->v, u->N * sizeof (double));
	u->w = (double *) realloc ((double *) u->w, u->N * sizeof (double));

	// the actual transformations into psi space
	cblas_dgemv (CblasRowMajor, CblasTrans, u->N, u->N,
			1.0, u->H, u->N, u->DT, 1, 0.0, u->v, 1);
	cblas_dgemv (CblasRowMajor, CblasTrans, u->N, u->N,
			1.0, u->H, u->N, u->ST, 1, 0.0, u->w, 1);
}

void time_step (hod * u)
{
	int i;

	for (i = 0; i <= u->N - 1; i++)
	{
		double v = u->v[i],
		       w = u->w[i],
		       a = u->E[i]*u->t;

		u->v[i] = v*cos (a) + w*sin (a);
		u->w[i] = w*cos (a) - v*sin (a);
	}

	// now we change them back into x-space
	cblas_dgemv (CblasRowMajor, CblasNoTrans, u->N, u->N,
			1.0, u->H, u->N, u->v, 1, 0.0, u->DT, 1);	
	cblas_dgemv (CblasRowMajor, CblasNoTrans, u->N, u->N,
			1.0, u->H, u->N, u->w, 1, 0.0, u->ST, 1);
}

void dump_step (hod * u)
{
	char * ani = (char *) malloc (30 * sizeof (char));

	// we will have to use a mkdir somewhere ...
	sprintf (ani, "%s/%05d.txt", u->dat, u->k);
	printf ("%s\n", ani);

	FILE * fout = fopen (ani, "w");
	free (ani);

	int i;
	for (i = 0; i <= u->N - 1; i++)
	{
		double x = u->h * (i - 0.5 * u->N),
		       P = pow(u->DT[i], 2) + pow(u->ST[i], 2),
		       V = pot (u, i);

		fprintf (fout, "% 15.8lf % 15.8lf % 15.8lf\n", x, P, V);
	}
	fprintf (fout, "\n");
	fclose (fout);
}

void create_frames (hod * u)
{
	// we create the appropriate directory for our
	// frames
	int status;
	status = mkdir (u->dat, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	dump_step (u);
	
	for (u->k = 1; u->k <= u->T-1; u->k++)
	{
		time_step (u);
		dump_step (u);
	}
}

void one_big_txt (hod * u)
{
	char * dat = (char *) malloc (20 * sizeof(char));

	sprintf (dat, "%s.txt", u->dat);

	FILE * fout = fopen (dat, "w");
	free (dat);

	int i;
	for (i = 0; i <= u->N - 1; i++)
	{
		double x = u->h * (i - 0.5 * u->N),
		       t = u->k * u->t;

		fprintf (fout, "% 15lf % 15.8lf % 15.8lf % 15.8lf\n",
				t, x, u->DT[i], u->ST[i]);
	}
	fprintf (fout, "\n");
	printf ("t = %d\n", u->k);

	for (u->k = 1; u->k <= u->T - 1; u->k++)
	{
		time_step (u);

		for (i = 0; i <= u->N - 1; i++)
		{
			double x = u->h * (i - 0.5 * u->N),
			       t = u->k * u->t;
			
			fprintf (fout, "% 15lf % 15.8lf % 15.8lf % 15.8lf\n",
					t, x, u->DT[i], u->ST[i]);
		}
		fprintf (fout, "\n");
		printf ("t = %d\n", u->k);
	}

	fclose (fout);
}

#endif

