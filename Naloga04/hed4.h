// hed4.h
////////////

#ifndef __HEADER_VRM4
#define __HEADER_VRM4

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct
{
	int N,               // vector dimension
	    T,               // temperature flag -- our time will be i/T
	    M,               // maximum time iteration
	    t,               // current integration index -- h*t = beta
	    s,               // which integrator to use ...
	    G;               // maximum number of vectors to be averaged

	gsl_rng * rand;      // random number generator

	FILE * fout;         // file for outputting the data

	double Z,            // vector of -z's for each psi we have a different value
	       F,            // free energy
	       h;            // step length -- temporal or thermal

	double complex ** g;  // vector of state vectors
} hod;

// funkcija zagrabi na i,i+1 mestu in ga transformira s
// propagatorjem u
void Utrans (hod * u, double complex a, int i, int j)
{
	int n = 2*i;

	double complex p  = a*u->h,
	               x1 = u->g[j][n+1];

	if (u->T == 1)
		p *= (-1)*I;

	// we take periodic boundaries 
	u->g[j][n]  *= cexp((-1)*p) * cexp(2*p);
	u->g[j][(n+3) % u->N] *= cexp ((-1)*p) * cexp (2*p);

	u->g[j][n+1] = cexp ((-1)*p) * (ccosh (2*p)*u->g[j][n+1] + csinh (2*p)*u->g[j][(n+2)%u->N]);
	u->g[j][(n+2)%u->N] = cexp ((-1)*p) * (ccosh (2*p)*u->g[j][(n+2)%u->N] + csinh (2*p)*x1);
}

void even (hod * u, double complex a, int j)
{
	int n = u->N/2,
	    i;

	for (i = 0; i <= n-1; i += 2)
		Utrans (u, a, i, j);
}

void odd (hod * u, double complex a, int j)
{
	int n = u->N/2,
	    i;

	for (i = 1; i <= n-1; i += 2)
		Utrans (u, a, i, j);
}

void Asym (hod * u)
{
	int j;
	for (j = 0; j <= u->G-1; j++)
	{
		even (u, 1.0 + 0.0*I, j);
		odd (u, 1.0 + 0.0*I, j);
	}
}

void S2doer (hod * u, double a, int j)
{
	even (u, a*0.5 + 0.0*I, j);
	odd (u, a*1.0 + 0.0*I, j);
	even (u, a*0.5 + 0.0*I, j);
}

void S3doer (hod * u, int con, double a, int j)
{
	double complex p1 = 0.125*(1 + sqrt(1.0/3)*I),
	       p5 = conj (p1),
	       p2 = 2*p1,
	       p4 = conj (p2),
	       p3 = 0.25 + 0.0*I;

	p1 *= a;
	p2 *= a;
	p3 *= a;
	p4 *= a;
	p5 *= a;

	if (con == 0)
	{
		even (u, p1, j);
		odd (u, p2, j);
		even (u, p3, j);
		odd (u, p4, j);
		even (u, p5, j);
	}

	else if (con == 1)
	{
		even (u, p5, j);
		odd (u, p4, j);
		even (u, p3, j);
		odd (u, p2, j);
		even (u, p1, j);
	}
}

void S2 (hod * u)
{
	int j;
	for (j = 0; j <= u->G-1; j++)
		S2doer (u, 1, j);
}

void S3 (hod * u)
{
	int j;
	for (j = 0; j <= u->G-1; j++)
		S3doer (u, 0, 1, j);
}

void S4 (hod * u)
{
	double f = pow (2, 1.0/3),
	       x0 = (-1)*f/(2 - f),
	       x1 = 1.0/(2 - f);

	int j;
	for (j = 0; j <= u->G-1; j++)
	{
		S2doer (u, x1, j);
		S2doer (u, x0, j);
		S2doer (u, x1, j);
	}
}

void S5 (hod * u)
{
	int j;
	for (j = 0; j <= u->G-1; j++)
	{
		S3doer (u, 0, 0.5, j);
		S3doer (u, 1, 0.5, j);
	}
}

void UpdateF (hod * u)
{
	u->F = (-1.0)/(u->t * u->h) * log (u->Z);
}

void UpdateZ (hod * u)
{
	u->Z = 0;

	int j;
	for (j = 0; j <= u->G-1; j++)
	{
		double S = 0;
		int i;
		for (i = 0; i <= u->N-1; i++)
			S += pow (cabs (u->g[j][i]),2);

		u->Z = u->Z * (1 - 1.0/(j+1)) + S/(j+1);
	}

	UpdateF (u);
}

void init_vec (hod * u, int j)
{
	u->g[j] = (double complex *) malloc (u->N * sizeof (double complex));

	int i;
	double S = 0;
	for (i = 0; i <= u->N-1; i++)
	{
		double x = gsl_ran_gaussian_ziggurat (u->rand, 2.0),
		       y = gsl_ran_gaussian_ziggurat (u->rand, 2.0);

		u->g[j][i] = x + y*I;
		S += pow (cabs(u->g[j][i]), 2);
	}

	S = sqrt(S);

	for (i = 0; i <= u->N-1; i++)
		u->g[j][i] /= S;
}

void init (hod * u, int N, int T, int M, int s, int G, double h, char * dat)
{
	u->N = pow(2,N);

	u->T = T;
	u->M = M;
	u->s = s;
	u->G = G;

	u->Z = 0;
	u->F = 0;
	u->h = h;

	u->fout = fopen (dat, "w");

	u->rand = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (u->rand, 10248023);

	int j;
	u->g = (double complex **) malloc (u->G * sizeof (double complex *));
	for (j = 0; j <= u->G-1; j++)
		init_vec (u, j);
}

void destroy (hod * u)
{
	int j;
	for (j = 0; j <= u->G-1; j++)
		free (u->g[j]);

	free (u->g);
	gsl_rng_free (u->rand);

	free (u);
}

void dumpZF (hod * u)
{
	fprintf (u->fout, "% 15lf % 15lf % 15lf\n", u->t*u->h, u->Z, u->F);
}

void simple_propagate (hod * u)
{

	switch (u->s)
	{
		case 1:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				Asym (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
		case 2:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				S2 (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
		case 3:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				S3 (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
		case 4:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				S4 (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
		case 5:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				S5 (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
		default:
			for (u->t = 0; u->t <= u->M-1; u->t++)
			{
				S4 (u);
				UpdateZ (u);
				dumpZF (u);
			}
			break;
	}
}

#endif
