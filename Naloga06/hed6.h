// hed6.h
////////////

#ifndef __HEADER_VRM6
#define __HEADER_VRM6

#include <stdbool.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

/* observables in 2D Ising's model */
typedef struct
{
	double M,			// <M>, magnetisation,
	       w,			// normalization of energy
	       E,			// <E>, energy
	       E2;			// <E^2>
} ising;

/* physical observables in Heisenberg's case */
typedef struct
{
	double * C,		// spin correlations: C(r) = <s_0 * s_r>
	       n;		// Heisenberg chain length
} chain;

/* the main struct ... it all depends on it */
typedef struct
{
	double * x;			// spin of the grid elements

	int n,				// grid dimensions, (n) in 1D and (n*n) in 2D case
	    v,				// annealing rate
	    I,				// current time index
	    J,				// (anti) fermomagnet parameter
	    max;			// maximum time iteration

	double h,			// magnetic field strength
	       T,			// temperature
	       dT,			// temperature step
	       dE,			// energy difference
	       H;			// E(sigma)

	FILE * fout,
	     * fani;
	
	char * basename;		// for plotting ...

	void * phys;			// physical observables -- either chain or ising

	gsl_rng * rand;			// random number generator

	bool (* step) (void *);
	void (* dump) (void *);
	void (* init_phys) (void *);
	void (* init_grid) (void *);
} hod;

void update_ising (hod * u, double sigma)
{
	double w1 = exp ((-1) * (u->H + u->dE)),
	       E = ((ising *) u->phys)->E,
	       E2 = ((ising *) u->phys)->E2,
	       w = ((ising *) u->phys)->w;

	u->H += u->dE;
	((ising *) u->phys)->M -= (2.0 * sigma)/(u->n * u->n);
	((ising *) u->phys)->E  = E/(1 + w) + ((u->H * w1)/(w + w1));
	((ising *) u->phys)->E2 = E2/(1 + w) + ((u->H * u->H * w1) / (w + w1));
	((ising *) u->phys)->w += w1;
}

bool step_ising (void * g)
{
	int n = ((hod *) g)->n,
	    J = ((hod *) g)->J,
	    h = ((hod *) g)->h,
	    i = (int) (n-1) * gsl_rng_uniform (((hod *) g)->rand),
	    j = (int) (n-1) * gsl_rng_uniform (((hod *) g)->rand);

	bool success = true;

	gsl_matrix_view m = gsl_matrix_view_array (((hod *) g)->x, n, n);

	double spin = gsl_matrix_get (&m.matrix, i, j),
	       s1   = gsl_matrix_get (&m.matrix, i, (j+1) % n),
	       s2   = gsl_matrix_get (&m.matrix, i, (j-1+n) % n),
	       s3   = gsl_matrix_get (&m.matrix, (i + 1) % n, j),
	       s4   = gsl_matrix_get (&m.matrix, (i - 1 + n) % n, j),
	       dE   = 2*spin * (h  + J*(s1 + s2 + s3 + s4));

	if (dE >= 0)
	{
	       double xi   = gsl_rng_uniform (((hod *) g)->rand);

	       if (xi > exp (((-1.0) * dE)/((hod *) g)->T))
		       success = false;
	       else
		       gsl_matrix_set (&m.matrix, i, j, (-1)*spin);
	}

	((hod *) g)->dE = dE;

	if (success)
		update_ising ((hod *) g, spin);

	return success;
}

void dump_ising_regular (hod * u)
{
	double chi = (1.0/u->T) * (1 - pow (((ising *) u->phys)->M, 2)),
	       C   = ((ising *) u->phys)->E2;

	C -= pow (((ising *) u->phys)->E, 2);
	C *= pow(1.0/u->T, 2);

	fprintf (u->fout, "% 12e % 12e % 12e\n", u->T, chi, C);
}

void dump_ising_animate (hod * u)
{
	int i, j;

	char * filename = (char *) malloc (40* sizeof (char));
	sprintf (filename, "%s/%05d.txt", u->basename, u->I);

	u->fani = fopen (filename, "w");

	for (i = 0; i <= u->n-1; i++)
	{
		for (j = 0; j <= u->n-1; j++)
			fprintf (u->fani, "% 2d", (int) u->x[j + i*u->n]);
		fprintf (u->fani, "\n");
	}

	fclose (u->fani);
	free (filename);
}

void dump_ising (void * u)
{
	dump_ising_regular ((hod *) u);
	dump_ising_animate ((hod *) u);
}

/* work in progress */
void update_chain (hod * u, double s1, double s2, double s3)
{
	return;
}

void dump_chain_regular (hod * u)
{
	return;
}

void dump_chain_animate (hod * u)
{
	return;
}

void dump_chain (void * u)
{
	dump_chain_regular ((hod *) u);
	dump_chain_animate ((hod *) u);
}

bool step_chain (void * u)
{
	/* until it's finished ... */
	return true;

	int i = (int) ((hod *) u)->n * gsl_rng_uniform (((hod *) u)->rand),
	    n = ((hod *) u)->n;

	bool success = true;

	// now we select the direction of rotation
	double costheta = (2 * gsl_rng_uniform (((hod *) u)->rand) - 1),
	       sintheta = sqrt (1 - costheta * costheta),
	       phi      = 2 * M_PI * gsl_rng_uniform (((hod *) u)->rand),
	       x1       = ((hod *) u)->x [i+1 + 0*n],
	       x2       = ((hod *) u)->x [i+1 + 1*n],
	       x3       = ((hod *) u)->x [i+1 + 2*n],
	       y1       = ((hod *) u)->x [i-1 + 0*n],
	       y2       = ((hod *) u)->x [i-1 + 1*n],
	       y3       = ((hod *) u)->x [i-1 + 2*n],
	       s1	= sin (phi) * sintheta,
	       s2	= sin (phi) * sintheta,
	       s3	= costheta,
	       J        = ((hod *) u)->J,
	       h        = ((hod *) u)->h,
	       dE       = (-1)*J*(s1*(x1 + y1) + s2*(x2 + y2) + s3*(x3 + y3)) - s3*h;

	if (dE >= 0)
	{
		double xi = gsl_rng_uniform (((hod *) u)->rand);

		if (xi > exp (((-1.0) * dE)/(((hod *) u)->T)))
			success = false;
		else
		{
			((hod *) u)->x [i+1 + 0*n] += s1;
			((hod *) u)->x [i+1 + 2*n] += s2;
			((hod *) u)->x [i+1 + 3*n] += s3;
		}
	}

	if (success)
		update_chain ((hod *) u, s1, s2, s3);

	return success;
}

void solver (hod * u)
{
	u->I = 0;
	u->T = 5;

	do
	{
		u->T -= u->dT;

		printf ("% 4d, % 10.3lf\n", u->I, u->T);

		int j = 0;
		do
		{
			bool success = u->step (u);

			if (success)
			{
				u->dump (u);
				j++;
			}
		} while (j <= u->v);

		u->I++;
	} while (u->T > 0);
}

/* let's just make them all point in one direction */
void init_grid_ising (void * g)
{
	int n = ((hod *) g)->n,
	    N2 = n*n,
	    i;

	((hod *) g)->x = (double *) malloc (N2 * sizeof (double));

	for (i = 0; i <= N2-1; i++)
		((hod *) g)->x [i] = 1;
}

void init_grid_chain (void * g)
{
	int n = ((hod *) g)->n,
	    N2 = 3 * n,
	    i;

	gsl_rng * rand = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rand, 5555555555);

	((hod *) g)->x = (double *) malloc (N2 * sizeof (double));

	for (i = 0; i <= n-1; i++)
	{
		// we get the random vector
		double x = 2 * gsl_rng_uniform (rand) - 1,
		       y = 2 * gsl_rng_uniform (rand) - 1,
		       z = 2 * gsl_rng_uniform (rand) - 1,
		       a = sqrt (x*x + y*y + z*z); // normalization

		// now that it's been normalized, we write it
		((hod *) g)->x [i+0*n] = x/a;
		((hod *) g)->x [i+1*n] = y/a;
		((hod *) g)->x [i+2*n] = z/a;
	}

	gsl_rng_free (rand);
}

void init_phys_ising (void * u)
{
	((hod *) u)->phys = (ising *) malloc (sizeof (ising));
	((ising *) (((hod *) u)->phys))->M  = 0;
	((ising *) (((hod *) u)->phys))->w  = 0;
	((ising *) (((hod *) u)->phys))->E  = 0;
	((ising *) (((hod *) u)->phys))->E2 = 0;

}

void init_phys_chain (void * u)
{
	((chain *) (((hod *) u)->phys))->n = ((hod *) u)->n;

	((chain *) (((hod *) u)->phys))->C =
		(double *) calloc (((hod *) u)->n, sizeof (double));
}

void init (hod * u,
		int mode,
		int n,
		int J,
		int v,
		int max,
		double dT,
		double h)
{
	u->n = n;
	u->J = J;
	u->v = v;
	u->max = max;

	u->T  = 0;
	u->dT = dT;
	u->h  = h;

	if (mode == 0)
	{
		u->step = &step_ising;
		u->dump = &dump_ising;
		u->init_grid = &init_grid_ising;
		u->init_phys = &init_phys_ising;
	}

	else
	{
		u->step = &step_chain;
		u->dump = &dump_chain;
		u->init_grid = &init_grid_chain;
		u->init_phys = &init_phys_chain;
	}

	u->basename = (char *) malloc (30 * sizeof (char));
	sprintf (u->basename, "c%d-n%d-J%d-h%d",
			mode, u->n, u->J, (int) u->h);

	mkdir (u->basename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	char * filename = (char *) malloc (34 * sizeof (char));
	sprintf (filename, "%s.txt", u->basename);

	u->fout = fopen (filename, "w");

	free (filename);

	u->rand = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (u->rand, 123123123);

	u->init_grid (u);
	u->init_phys (u);
}

void destroy (hod * u)
{
	free (u->x);
	free (u->phys);
	free (u->basename);
	gsl_rng_free (u->rand);
	fclose (u->fout);

	free (u);
}

#endif

