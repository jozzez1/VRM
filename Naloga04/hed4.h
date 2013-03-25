// hed4.h
////////////

#ifndef __HEADER_VRM4
#define __HEADER_VRM4

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// our "matrix" for matrix element calculation
typedef struct
{
	int ** c,       // connections array
	    * b,        // number in binary
	    N,          // total number of bits
	    m,          // number of ones
	    n;          // number of connections

	// c[0][j] -- connection factor,
	// c[1][j] -- connects to, for j-th connection
} pov;

// the walker struct
typedef struct
{
	int N,                // vector dimension
	    T,                // temperature flag -- our time will be i/T
	    M,                // maximum time iteration
	    t,                // current integration index -- h*t = beta
	    c,                // calculate-hamiltonian-flag
	    E,                // yes to energy!
	    C,                // time correlation of the first spin
	    J,                // yes to spin current! :PPPP
	    G;                // maximum number of vectors to be averaged

	gsl_rng * rand;       // random number generator

	FILE * fout;          // file for outputting the data

	void (* connect) (pov *, int); // function pointer to connection function
	void (* stepper) (void *);     // function pointer propagator

	double Z,             // vector of -z's for each psi we have a different value
	       F,             // free energy
	       d,             // anisotropic factor in the hamiltonian
	       h;             // step length -- temporal or thermal

	pov ** e;             // vector of state vecor connections

	double complex ** g,  // vector of state vectors
	       H;             // Hamiltonian average energy/Spin current/1st spin in z
} hod;

// propagator U grabs on the i,i+1 component of the j-th vector
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

void Asym (void * u)
{
	int j;
	for (j = 0; j <= ((hod *)u)->G-1; j++)
	{
		even ((hod *) u, 1.0 + 0.0*I, j);
		odd ((hod *) u, 1.0 + 0.0*I, j);
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

void S2 (void * u)
{
	int j;
	for (j = 0; j <= ((hod *) u)->G-1; j++)
		S2doer ((hod *) u, 1, j);
}

void S3 (void * u)
{
	int j;
	for (j = 0; j <= ((hod *) u)->G-1; j++)
		S3doer ((hod *) u, 0, 2, j); // used to be u,0,1,j
}

void S4 (void * u)
{
	double f = pow (2, 1.0/3),
	       x0 = (-1)*f/(2 - f),
	       x1 = 1.0/(2 - f);

	int j;
	for (j = 0; j <= ((hod *) u)->G-1; j++)
	{
		S2doer ((hod *) u, x1, j);
		S2doer ((hod *) u, x0, j);
		S2doer ((hod *) u, x1, j);
	}
}

void S5 (void * u)
{
	int j;
	for (j = 0; j <= ((hod *) u)->G-1; j++)
	{
		                     // fixed for correct time
		S3doer ((hod *) u, 0, 1, j); // used to be u,0,0.5,j
		S3doer ((hod *) u, 1, 1, j); // used to be u,1,0.5,j
	}
}

void UpdateF (hod * u)
{
	u->F = (-2.0)/(u->t * u->h) * log (u->Z);
}

// takes care of diagonal connections
void binary (pov * e, int x)
{
	int i;
	e->c[0][0] = 0;
	e->c[1][0] = x;
	for (i = 0; i <= e->N-1; i++)
	{
		e->b[i] = x%2;
		x /= 2;

		if (e->b[i] == 1)
			e->m++;
	}
}

// identity connection
void connect_Id (pov * e, int x)
{
	e->c[0][0] = 1;
}

// matrix elements of Hamiltonian
void connect_E (pov * e, int x)
{
	int i;
	for (i = 0; i <= e->N - 1; i++)
	{
		// if we find 01, but b is reversed, thus 10
		if (e->b[i] == 1 && e->b[(i+1)%e->N] == 0)
		{
			e->n++;
			e->c[0] = (int *) realloc ((int *) e->c[0],
					e->n * sizeof (int));
			e->c[1] = (int *) realloc ((int *) e->c[1],
					e->n * sizeof (int));

			e->c[0][e->n-1] = 2;
			e->c[1][e->n-1] = x - (int) pow (2, i) + (int) pow (2, (i+1)%e->N);

			e->c[0][0] -= 1;
		}

		// if we find an 10, but b is reversed, thus 01
		if (e->b[i] == 0 && e->b[(i+1)%e->N] == 1)
		{
			e->n++;
			e->c[0] = (int *) realloc ((int *) e->c[0],
					e->n * sizeof (int));
			e->c[1] = (int *) realloc ((int *) e->c[1],
					e->n * sizeof (int));

			e->c[0][e->n-1] = 2;
			e->c[1][e->n-1] = x - (int) pow (2, (i+1)%e->N) + (int) pow (2, i);

			e->c[0][0] -= 1;
		}

		if ((e->b[i] == 0 && e->b[(i+1)%e->N] == 0) || 
				(e->b[i] == 1 && e->b[(i+1)%e->N] == 1))
			e->c[0][0] += 1;
	}
}

// matrix elements of J -- the spin currrent
void connect_J (pov * e, int x)
{
	int i;
	for (i = 0; i <= e->N - 1; i++)
	{
		if (e->b[i] == 1 && e->b[(i+1)%e->N] == 0)
		{
			e->n++;
			e->c[0] = (int *) realloc ((int *) e->c[0],
					e->n * sizeof (int));
			e->c[1] = (int *) realloc ((int *) e->c[1],
					e->n * sizeof (int));

			e->c[0][e->n-1] = -2;
			e->c[1][e->n-1] = x - (int) pow (2, i) + (int) pow (2, (i+1)%e->N);
		}

		if (e->b[i] == 0 && e->b[(i+1)%e->N] == 1)
		{
			e->n++;
			e->c[0] = (int *) realloc ((int *) e->c[0],
					e->n * sizeof (int));
			e->c[1] = (int *) realloc ((int *) e->c[1],
					e->n * sizeof (int));

			e->c[0][e->n-1] = 2;
			e->c[1][e->n-1] = x - (int) pow (2, (i+1)%e->N) + (int) pow (2, i);
		}
	}
}

// matrix elements of s^{z}_{1}
void connect_C (pov * e, int x)
{
	e->c[0][0] = e->b[e->N-1];
	if (!e->c[0][0])
		e->c[0][0]--;
}

// init e
void init_e (pov * e, int N)
{
	e->N = N;
	e->m = 0;
	e->n = 1;

	e->b = (int *) malloc (e->N * sizeof (int));

	e->c = (int **) malloc (2 * sizeof (int *));
	
	e->c[0] = (int *) malloc (sizeof (int));
	e->c[1] = (int *) malloc (sizeof (int));
}

void destroy_e (pov * e)
{
	free (e->b);
	free (e->c[0]);
	free (e->c[1]);
	free (e->c);
	free (e);
}

void UpdateH (hod * u)
{
	u->H = 0.0 + 0.0 * I;

	int i, j, k;
	for (j = 0; j <= u->G-1; j++)
	{
		double complex S = 0.0 + 0.0*I;
		for (i = 0;  i <= u->N-1; i++)
		{
			for (k = 0; k <= u->e[i]->n-1; k++)
				S += conj (u->g[j][i]) * u->e[i]->c[0][k] * u->g[j+u->T*u->G][u->e[i]->c[1][k]];

			if (u->T)
				S *= cpow (I, u->e[i]->n-1);
		}

		u->H = u->H * (1 - 1.0/(j+1)) + S/(j+1);
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

	if (u->E || u->C || u->J)
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

void spin_corr_init_vec (hod * u)
{
	int i, j, k,
	    g = u->G;

	// we extend the |psi> for X|psi> = |chi>
	u->G *= 2;
	u->g = (double complex **) realloc ((double complex **) u->g, u->G * sizeof (double complex *));

	// we initialize the |chi>-s
	for (j = g; j <= u->G-1; j++)
	{
		u->g[j] = (double complex *) malloc (u->N * sizeof (double complex));
		for (i = 0; i <= u->N - 1; i++)
		{
			for (k = 0; k <= u->e[i]->n - 1; k++)
				u->g[j][i] = u->g[j-g][u->e[i]->c[1][k]] * u->e[i]->c[0][k];
		}
	}
	u->G = g;
}

// we "plot" our connectivity as a matrix
void controlH (hod * u)
{
	int N = u->N,
	    i, j;

	int ** H = (int **) malloc (N * sizeof(int *));
	for (i = 0; i <= N-1; i++)
	{
		H[i] = (int *) calloc (N, sizeof(int));
		for (j = 0; j <= u->e[i]->n-1; j++)
			H[i][u->e[i]->c[1][j]] += u->e[i]->c[0][j];

		if (u->e[i]->N <= 5)
		{
			printf ("|");
			for (j = 0; j <= N-1; j++)
			{
				if (H[i][j] != 0)
					printf ("% 2d", H[i][j]);
				else
					printf ("  ");
			}		
			printf ("|\n");
		}
		
	}

	for (i = 0; i <= N-1; i++)
		free (H[i]);
	free (H);

	exit (EXIT_SUCCESS);
}

void init (hod * u, int N, int T, int E, int C, int J, int M, int s, int G, int c, double h, char * dat)
{
	u->N = pow(2,N);

	u->E = E;
	u->C = C;
	u->J = J;
	u->T = T;
	u->M = M;
	u->G = G;
	u->c = c;

	u->Z = 0;
	u->F = 0;
	u->h = h;

	u->H = 0.0 + 0.0*I;

	u->fout = fopen (dat, "w");

	u->rand = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (u->rand, 10248023);

	int j;
	u->g = (double complex **) malloc (u->G * sizeof (double complex *));
	for (j = 0; j <= u->G-1; j++)
		init_vec (u, j);

	if (u->E) u->connect = &connect_E;
	else if (u->C) u->connect = &connect_C;
	else if (u->J) u->connect = &connect_J;
	else u->connect = &connect_Id;

	// we select the stepper
	switch (s)
	{
		case 1:
			u->stepper = &Asym;
			break;
		case 2:
			u->stepper = &S2;
			break;
		case 3:
			u->stepper = &S3;
			break;
		case 4:
			u->stepper = &S4;
			break;
		case 5:
			u->stepper = &S5;
			break;
		default:
			u->stepper = &S4;
			break;
	}
	
	u->e = (pov **) malloc (u->N * sizeof (pov *));
	for (j = 0; j <= u->N-1; j++)
	{
		u->e[j] = (pov *) malloc (sizeof (pov));
		init_e (u->e[j], N);
		binary (u->e[j], j);
		u->connect (u->e[j], j);
	}

	if (u->T) spin_corr_init_vec (u);
	if (u->c)
		controlH (u);
}

void destroy (hod * u)
{
	int j;
	for (j = 0; j <= u->G*(1 + u->T)-1; j++)
		free (u->g[j]);

	if (u->E)
	{
		for (j = 0; j <= u->N-1; j++)
			destroy_e (u->e[j]);

		free (u->e);
	}

	free (u->g);
	gsl_rng_free (u->rand);
	
	fclose (u->fout);

	free (u);
}

void dump (hod * u)
{
	if (!u->E && !u->C && !u->J)
		fprintf (u->fout, "% 15lf % 15e % 15e\n", u->t*u->h, u->Z, u->F);

	else if (u->E || u->C || u->J)
		fprintf (u->fout, "% 15lf % 15e % 15e % 15e % 15e\n",
				u->t*u->h, u->Z, u->F, creal (u->H), cimag (u->H));

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

void simple_propagate (hod * u)
{
	for (u->t = 0; u->t <= u->M-1; u->t++)
	{
		u->stepper (u);
		UpdateZFH (u);
		dump (u);
	}
}

#endif

