// hed1.h
////////////

#ifndef __HEADER_VRM1
#define __HEADER_VRM1

#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

typedef struct
{
	gsl_vector_complex * psi0, // old iteration
			   * psi1, // next iteration
			   * bpr;  // between result

	gsl_matrix_complex * A,    // the LU decomposition, and matri A
			   * B;    // the adjoint matrix of A
	
	gsl_permutation * p;       // permuation for LU decomposition

	FILE * fout;               // pointer to a file

	char * dat;                // "animation" directory

	double h,                  // space step
	       t,                  // time step
	       a,                  // shift from the equilibrium
	       Lambda;             // parameter in the potential

	int N,                     // maximum time step
	    M,                     // space granularity
	    pet,                   // flag for 5-diagonal matrix
	    n;                     // current time iteration

} hod;

// a.k.a. data dumer
void dump (hod * u)
{
	int i;
	for (i = 0; i <= u->M-1; i++)
	{
		gsl_complex a = gsl_vector_complex_get (u->psi0, i);
		fprintf (u->fout, "% 15lf % 15lf % 15lf % 15lf\n",
				u->n*u->t, u->h*(i - u->M/2), a.dat[0], a.dat[1]);
	}
	fprintf (u->fout, "\n");
}

// shit out a portion of data for animation
void shit (hod * u)
{
	char * dat1 = (char *) malloc (65 * sizeof(char));
	sprintf (dat1, "%s/%05d.txt", u->dat, u->n);

	printf ("%s\n", dat1);
	u->fout = fopen (dat1, "w");

	int i;
	for (i = 0; i <= u->M-1; i++)
	{
		gsl_complex a = gsl_vector_complex_get (u->psi0, i);

		double r = a.dat[0]*a.dat[0] + a.dat[1]*a.dat[1],
		       x = u->h* (i - u->M/2),
		       V = 0.5*x*x + u->Lambda * pow (x, 4);

		fprintf (u->fout, "% 15lf % 15lf % 15lf\n", u->h*(i - u->M/2), r, V);
	}
	fprintf (u->fout, "\n");

	fclose (u->fout);
}

// potential calc
double vofm (int m, hod * u)
{
	double x = m * u->h - u->M * u->h/2, // we shift the potential by 1/2
	       L = u->Lambda,
	       V = 0.5*x*x + L * pow (x, 4);

	return V;
}

// we prepare the starting vectors
// N=0
void start0 (hod * u)
{
	int n;

	gsl_complex phi;
	phi.dat[1] = 0.0;

	for (n = 1; n <= u->M-2; n++)
	{
		double x = u->h*(n - u->M/2) - u->a;

		phi.dat[0] = exp((-1)*x*x/2)/pow(M_PI, 0.25);

		gsl_vector_complex_set (u->psi0, n, phi);
	}

	phi.dat[0] = 0.0;

	gsl_vector_complex_set (u->psi0, 0, phi);
	gsl_vector_complex_set (u->psi0, u->M-1, phi);
}

// N=1
void start1 (hod * u)
{
	int n;

	gsl_complex phi;
	phi.dat[1] = 0.0;

	for (n = 1; n <= u->M-2; n++)
	{
		double x = u->h*(n - u->M/2) - u->a;

		phi.dat[0] = 2*x*exp((-1)*x*x/2)/(sqrt(2)*pow(M_PI, 0.25));
		
		gsl_vector_complex_set (u->psi0, n, phi);
	}

	phi.dat[0] = 0.0;

	gsl_vector_complex_set (u->psi0, 0, phi);
	gsl_vector_complex_set (u->psi0, u->M-1, phi);
}

// N=2
void start2 (hod * u)
{
	int n;

	gsl_complex phi;
	phi.dat[1] = 0.0;

	for (n = 1; n <= u->M-2; n++)
	{
		double x = u->h*(n - u->M/2) - u->a;

		phi.dat[0] = (4*x*x - 2)*exp((-1)*x*x/2)/(sqrt(8)*pow(M_PI, 0.25));
		
		gsl_vector_complex_set (u->psi0, n, phi);
	}

	phi.dat[0] = 0.0;

	gsl_vector_complex_set (u->psi0, 0, phi);
	gsl_vector_complex_set (u->psi0, u->M-1, phi);
}

void matrixA (hod * u)
{
	int i;

	gsl_complex b;
	b.dat[1] = (-1)*u->t*0.25/(u->h * u->h);
	b.dat[0] = 0;

	// the center ...
	for (i = 1; i <= u->M-2; i++)
	{
		double V = vofm (i, u);

		gsl_complex a;

		a.dat[0] = 1;
		a.dat[1] = u->t/(2 * u->h * u->h)*(1 + u->h * u->h * V);

		gsl_matrix_complex_set (u->A, i, i, a);
		gsl_matrix_complex_set (u->A, i, i+1, b);
		gsl_matrix_complex_set (u->A, i, i-1, b);
	}

	// now the margins
	gsl_matrix_complex_set (u->A, 0, 1, b);
	gsl_matrix_complex_set (u->A, u->M-1, u->M-2, b);
	
	b.dat[1] = u->t/(u->h * u->h)*(2 + u->h * u->h * vofm(0, u));
	gsl_matrix_complex_set (u->A, 0, 0, b);

	b.dat[1] = u->t/(u->h * u->h)*(2 + u->h * u->h * vofm(u->M-1, u));
	gsl_matrix_complex_set (u->A, u->M-1, u->M-1, b);

}

void matrixA1 (hod * u)
{
	int i;

	gsl_complex b1 = gsl_complex_rect (0, (-1)*(u->t/(3.0*u->h*u->h)));
	gsl_complex b2 = gsl_complex_rect (0, u->t/(u->h*u->h*48));

	// center ...
	for (i = 2; i <= u->M-3; i++)
	{
		double V = vofm (i, u);

		gsl_complex a;

		a.dat[0] = 1;
		a.dat[1] = (u->t/(4*u->h*u->h))*(2.5 + 2.0*u->h*u->h*V);

		gsl_matrix_complex_set (u->A, i, i-1, b1);
		gsl_matrix_complex_set (u->A, i, i-2, b2);
		gsl_matrix_complex_set (u->A, i, i, a);
		gsl_matrix_complex_set (u->A, i, i+1, b1);
		gsl_matrix_complex_set (u->A, i, i+2, b2);
	}

	// and now margins ...
	gsl_matrix_complex_set (u->A, 0, 1, b1);
	gsl_matrix_complex_set (u->A, 0, 2, b2);
	gsl_matrix_complex_set (u->A, 1, 0, b1);
	gsl_matrix_complex_set (u->A, 1, 2, b1);
	gsl_matrix_complex_set (u->A, 1, 3, b2);

	gsl_complex a;
	a.dat[0] = 1;
	a.dat[1] = (u->t/(4*u->h*u->h))*(2.5 + 2.0*u->h*u->h*vofm(0, u));
	gsl_matrix_complex_set (u->A, 0, 0, a);
	a.dat[1] = (u->t/(4*u->h*u->h))*(2.5 + 2.0*u->h*u->h*vofm(1, u));
	gsl_matrix_complex_set (u->A, 1, 1, a);
	a.dat[1] = (u->t/(4*u->h*u->h))*(2.5 + 2.0*u->h*u->h*vofm(u->M-2, u));
	gsl_matrix_complex_set (u->A, u->M-2, u->M-2, a);
	a.dat[1] = (u->t/(4*u->h*u->h))*(2.5 + 2.0*u->h*u->h*vofm(u->M-1, u));
	gsl_matrix_complex_set (u->A, u->M-2, u->M-1, a);

	gsl_matrix_complex_set (u->A, u->M-1, u->M-2, b1);
	gsl_matrix_complex_set (u->A, u->M-1, u->M-3, b2);
	gsl_matrix_complex_set (u->A, u->M-2, u->M-1, b1);
	gsl_matrix_complex_set (u->A, u->M-2, u->M-3, b1);
	gsl_matrix_complex_set (u->A, u->M-2, u->M-4, b2);
}

// construct B from A
void herm_adA (hod * u)
{
	int i;
	gsl_complex a;

	for (i = 0; i <= u->M-1; i++)
	{
		a = gsl_complex_conjugate(gsl_matrix_complex_get (u->A, i, i));
		gsl_matrix_complex_set (u->B, i, i, a);
	}

	a = gsl_complex_conjugate (gsl_matrix_complex_get (u->A, 0, 1));

	for (i = 0; i <= u->M-2; i++)
		gsl_matrix_complex_set (u->B, i, i+1, a);
	for (i = 1; i <= u->M-1; i++)
		gsl_matrix_complex_set (u->B, i, i-1, a);

	// in case we have additional diagonals
	if (u->pet == 1)
	{
		a = gsl_complex_conjugate (gsl_matrix_complex_get (u->A, 0, 2));
		for (i = 0; i <= u->M-3; i++)
			gsl_matrix_complex_set (u->B, i, i+2, a);
		for (i = 2; i <= u->M-1; i++)
			gsl_matrix_complex_set (u->B, i, i-2, a);
	}
}

// initialization function
void init (hod * u, double h, double t, double a, double Lambda,
		int N, int M, int switsch, int pet)
{
	u->N = N;
	u->M = M;
	u->n = 0;
	u->pet = pet;

	u->t = t;
	u->a = a;
	u->Lambda = Lambda;
	u->h = h;

	u->dat = (char *) malloc (40 * sizeof(char));

	// we allocate the space
	u->p = gsl_permutation_alloc (u->M);

	u->A = gsl_matrix_complex_calloc (u->M, u->M);
	u->B = gsl_matrix_complex_calloc (u->M, u->M);

	u->psi0 = gsl_vector_complex_calloc (u->M);
	u->psi1 = gsl_vector_complex_calloc (u->M);
	u->bpr = gsl_vector_complex_calloc (u->M);

	// we initialize the matrices
	// initialize A
	if (u->pet == 0)
		matrixA (u);
	if (u->pet == 1)
		matrixA1 (u);

	herm_adA (u);     // initialize B

	// we decompose the matrix A
	int signum; // unnecessary for solving the problem
	gsl_linalg_complex_LU_decomp (u->A, u->p, &signum); // decomposition

	// starting vector initialization
	switch (switsch)
	{
		case 0:
			start0 (u);
			break;
		case 1:
			start1 (u);
			break;
		case 2:
			start2 (u);
			break;
		default:
			printf ("Choosing ground state\n");
			start0 (u);
			break;
	}
}

void swap (gsl_vector_complex * x, gsl_vector_complex * y)
{
	int N = x->size;
	gsl_vector_complex * z = gsl_vector_complex_alloc (N);

	int i;
	for (i = 0; i <= N-1; i++)
	{
		gsl_vector_complex_set (z, i, gsl_vector_complex_get(x, i));
		gsl_vector_complex_set (x, i, gsl_vector_complex_get(y, i));
		gsl_vector_complex_set (y, i, gsl_vector_complex_get(z, i));
	}

	gsl_vector_complex_free (z);
}

void step (hod * u)
{
	gsl_complex a = gsl_complex_rect (1.0, 0.0);
	gsl_complex b = gsl_complex_rect (0.0, 0.0);

	gsl_blas_zgemv (CblasNoTrans, a, u->B, u->psi0, b, u->bpr);
	gsl_linalg_complex_LU_solve (u->A, u->p, u->bpr, u->psi1);
}

void solver (hod * u)
{
	dump (u);

	for (u->n = 1; u->n <= u->N-1; u->n++)
	{
		gsl_complex a = gsl_complex_rect (0, 0);

		step (u);
		swap (u->psi0, u->psi1);

		gsl_vector_complex_set (u->psi0, 0, a);
		gsl_vector_complex_set (u->psi0, u->M-1, a);

		dump (u);
	}
}

void shitter (hod * u)
{
	shit (u);
	for (u->n = 1; u->n <= u->N-1; u->n++)
	{
		gsl_complex a = gsl_complex_rect (0, 0);

		step (u);
		swap (u->psi0, u->psi1);

		gsl_vector_complex_set (u->psi0, 0, a);
		gsl_vector_complex_set (u->psi0, u->M-1, a);

		shit (u);
	}
}

void destroy (hod * u)
{
	gsl_vector_complex_free (u->psi0);
	gsl_vector_complex_free (u->psi1);

	gsl_matrix_complex_free (u->A);
	gsl_matrix_complex_free (u->B);

	gsl_permutation_free (u->p);

	free (u);
}

#endif
