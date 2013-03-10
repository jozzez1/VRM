// hed2.h
////////////

#ifndef __HEADER_VRM2
#define __HEADER_VRM2

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <cblas.h>
#include <clapack.h>

typedef struct
{
	double h,
	       L;

	int N,
	    p;       // red odvodov

	double * H,  // Hamiltonian, oz. na koncu matrika v
	       * v,  // eigenvector matrix in Lanczos space
	       * w,  // Lanczos transformation matrix
	       * DT, // dominant of tridiagonal
	       * ST, // subdiagonal of tridiagonal
	       * e;  // eigenvalue array
} hod;

// potential
double pot (hod * u, int i)
{
	double x = u->h * (i - 0.5*u->M),
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

	gsl_matrix view m = gsl_matrix_view_array (u->H, u->N, u->N);
gg
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

void swapv (double * v1, double * v2, int N)
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

	gsl_matrix_view a = gsl_matrix_view_array (u->w, u->N, u->N);

	double * b = (double *) malloc (u->N * sizeof (double)),
	       * a = (double *) malloc (u->N * sizeof (double)),
	       * v0 = (double *) calloc (u->N, sizeof (double)),
	       * v1 = (double *) calloc (u->N, sizeof (double)),
	       * w = (double *) calloc (u->N, sizeof (double));

	b[0] = 0;
	v1[0] = 1;

	for (i = 0; i <= u->N - 2; i++)
		gsl_matrix_set (&a.matrix, i, 0, v1[i]);

	for (i = 0; i <= u->N-2; i++)
	{
		cblas_dsymv (CblasRowMajor, CblasUpper, u->N, 1.0, u->H, u->N, v1, 1, 0.0, w, 1);
		a[i] = scalar (w, v1, u->N);

		int j;
		for (j = 0; j <= u->N - 1; j++)
			w[j] -= (a[i]*v1[j] + b[i]*v0[j]);

		b[i+1] = sqrt (scalar (w, w, u->N));

		swap (v1, v0);

		for (j = 0; j <= u->N-1; j++)
		{
			v1[j] = w[j]/b[i+1];
			gsl_matrix_set (&a.matrix, j, i, v1[j]);
		}
	}

	cblas_dsymv (CblasRowMajor, CblasUpper, u->N, 1.0, u->H, u->N, v1, 1, 0.0, w, 1);
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
void init (hod * u, double h, double L, int N, int p)
{
	u->h = h;
	u->L = L;
	u->N = N;
	u->p = p;

	int M = u->N * u->N,
	    D = u->N - 1;

	u->H = (double *) calloc (M, sizeof(double));
	u->v = (double *) malloc (M * sizeof (double));
	u->DT = (double *) malloc (u->N * sizeof (double));
	u->ST = (double *) malloc (D * sizeof (double));
	u->e = (double *) malloc (u->N * sizeof (double));

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
	free (H);
	free (v);
	free (DT);
	free (ST);
	free (e);

	free (u);
}


// now to diagonalize this baby with this handy wrapper ...
// compile with: -llapack -lblas -lcblas -latlas
// for fortran code ALL arguments must be passed as pointers
static long MRRR (hod * u)
{
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
			int * tryrac,     // (in/out) logical (int or long int) ... for high accuracy
			double * work,    // (out) array of dimension lwork, whic is the output
			int * lwork,      // (in) for JOBZ='V' has to be >= max (1, 12*N)
			double * iwork,   // (out) same as work
			double * liwork,  // (in) for JOBZ='V' has to be >= max (1, 10*N)
			int * Info);      // (out) =0 for successfull operation

	char JOBZ   =   'V',
	     RANGE  =   'A';

	int i = u->N-1;

	double * d = (double *) malloc (u->N * sizeof(double)),
	       * e = (double *) malloc (i * sizeof (double));
	
	for (i = 0; i <= u->N - 1; i++)
	{
		d[i] = u->DT[i];
		e[i] = u->ST[i];
	}

	double VL = 0,
	       VU = 1;

	int IL = 0,
	    IU = 1,
	    M;

	// for 'w' I will use 'u->e' -- already an array
	// for 'z' I will use 'u->v' -- already an array
	
	int LDZ    = u->N,
	    NZC    = u->N,
	    ISUPPZ = (int *) malloc (2*u->N * sizeof (int)),
	    TRYRAY = 1, // I hope this means 'yes, check for precision :D'
	    LWORK  = 12*u->N,
	    LIWORK = 10*u->N,
	    INFO;

	double * WORK  = (double *) malloc (LWORK * sizeof (double)),
	       * IWORK = (double *) malloc (LIWORK * sizeof (double));


	dstemr_ (&JOBZ, &RANGE, &u->N, d, e, &VL, &VU, &IL, &IU, &M, u->e, u->v, &LDZ, &NZC,
			ISUPPZ, &TRYRAY, WORK, &LWORK, IWORK, &LIWORK, &INFO);

	free (d);
	free (e);
	free (ISUPPZ);
	free (WORK);
	free (IWORK);

	return INFO;
}

#endif
