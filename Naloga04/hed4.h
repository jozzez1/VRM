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
	    E,               // yes to energy!
	    G;               // maximum number of vectors to be averaged

	gsl_rng * rand;      // random number generator

	FILE * fout;         // file for outputting the data

	double Z,            // vector of -z's for each psi we have a different value
	       F,            // free energy
	       h;            // step length -- temporal or thermal

	double complex ** g,  // vector of state vectors
	       H,             // Hamiltonian average energy -- can be a complex number
	       C;             // Spin correlation over time
} hod;

// propagator U grabs on the i,i+1 component of the j-th
// vector
void Utrans (hod * u, double complex a, int i, int j)
{
	int n = 2*i;

	double complex p  = (-1)*a*u->h,
	               x1 = u->g[j][n+1];

	// if we will use time
	if (u->T == 2)
		p *= I;

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
	u->F = (-2.0)/(u->t * u->h) * log (u->Z);
}

// with this we do <i,i+1| H |i,i+1>_j ... to make it clear
// it's highly optimized
void UpdateH (hod * u)
{
	u->H = 0.0 + 0.0*I;

	int n = u->N/2,
	    i, j;
	for (j = 0; j <= u->G-1; j++)
	{
		double complex S = 0.0 + 0.0*I;
		for (i = 1; i <= n; i++)
			S += conj (u->g[j][(2*i)%u->N]) * u->g[j][(2*i-1)];

		S *= 2;
		u->H = u->H*(1 - 1.0/(j+1)) + S/(j+1);
	}

	u->H /= u->Z;
}

void UpdateZFH (hod * u)
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

	if (u->E == 1)
		UpdateH (u);
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

void init (hod * u, int N, int T, int E, int M, int s, int G, double h, char * dat)
{
	u->N = pow(2,N);

	u->E = E;
	u->T = T;
	u->M = M;
	u->s = s;
	u->G = G;

	u->Z = 0;
	u->F = 0;
	u->h = h;

	u->H = 0.0 + 0.0*I;
	u->C = 0.0 + 0.0*I;

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
	
	fclose (u->fout);

	free (u);
}

void dump (hod * u)
{
	if (u->E == 0 && u->T == 0)
		fprintf (u->fout, "% 15lf % 15e % 15e\n", u->t*u->h, u->Z, u->F);

	else if (u->E == 1 && u->T == 0)
		fprintf (u->fout, "% 15lf % 15e % 15e % 15e % 15e\n",
				u->t*u->h, u->Z, u->F, creal (u->H), cimag (u->H));

	else if (u->T == 2)
		fprintf (u->fout, "% 15lf % 15e % 15e\n",
				u->t*u->h, creal (u->C), cimag (u->H));

	// progress bar
	int percent = 20*u->t/u->M,
	    i;

	printf ("% 4.0lf%% [", (100.0*u->t)/u->M);
	for (i = 0; i <= percent-1; i++)
		printf ("=");
	printf (">");
	for (; i<= 18; i++)
		printf (" ");
	printf ("]\n\033[F\033[J");
}

/* time correlations */
void spin_corr_init_vec (hod * u)
{
	int g = u->G,
	    n = u->N/2,
	    i, j;

	// we extend the |psi> for sigma^z |psi> = |chi>
	u->G *= 2;
	u->g = (double complex **) realloc ((double complex **) u->g, u->G * sizeof (double complex *));

	// we initialize the |chi>-s
	for (j = g; j <= u->G-1; j++)
	{
		u->g[j] = (double complex *) malloc (u->N * sizeof(double complex));
		for (i = 0; i <= n-1; i++)
		{
			u->g[j][2*i] = u->g[j-g][2*i];
			u->g[j][2*i+1] = (-1)*u->g[j-g][2*i+1];
		}
	}
}

/* spin corr update function */
void UpdateC (hod * u)
{
	u->C = 0.0 + 0.0*I;

	int n = u->N/2,
	    g = u->G/2,
	    i, j;

	for (j = 0; j <= g-1; j++)
	{
		double complex S = 0.0 + I*0.0;
		for (i = 0; i <= n-1; i++)
		{
			S += conj (u->g[j][2*i]) * u->g[j+g][2*i];
			S -= conj (u->g[j][2*i + 1]) * u->g[j+g][2*i + 1];
		}

		u->C = u->C * (1 - 1.0/(j+1)) + S/(j+1);
	}
}

void simple_propagate (hod * u)
{
	if (u->T == 0)
	{
		switch (u->s)
		{
			case 1:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					Asym (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
			case 2:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S2 (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
			case 3:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S3 (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
			case 4:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S4 (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
			case 5:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S5 (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
			default:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S4 (u);
					UpdateZFH (u);
					dump (u);
				}
				break;
		}
	}

	else if (u->T == 2)
	{
		spin_corr_init_vec (u);
		switch (u->s)
		{
			case 1:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					Asym (u);
					UpdateC (u);
					dump (u);
				}
				break;
			case 2:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S2 (u);
					UpdateC (u);
					dump (u);
				}
				break;
			case 3:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S3 (u);
					UpdateC (u);
					dump (u);
				}
				break;
			case 4:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S4 (u);
					UpdateC (u);
					dump (u);
				}
				break;
			case 5:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S5 (u);
					UpdateC (u);
					dump (u);
				}
				break;
			default:
				for (u->t = 0; u->t <= u->M-1; u->t++)
				{
					S4 (u);
					UpdateC (u);
					dump (u);
				}
				break;
		}
	}
}

#endif

