// hed6.h
////////////

#ifndef __HEADER_VRM6
#define __HEADER_VRM6

#include <sys/stat.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>


// threading for plotter
#include <omp.h> // !!! requires the GNU Compiler Collection (gcc) !!! -- already in the make file

static const double TC = 2.269185314;

/* we start over ... */
typedef struct
{
	int n,		// grid size = n*n
	    J,		// well some parameter, you can figure out which
	    I,		// internal index
	    v,		// averaging time
	    mode,	// if we measure hysteresis (after heating, we cool the system)
	    k;		// internal index

	double h,	// magnetic field strength
	       H,	// current energy
	       dE,	// energy difference
	       E,	// <E>
	       E2,	// <E^2>
	       M,	// <M>
	       m,	// average spin <S>
	       T,	// 1/beta
	       dT,	// temperature step
	       Tmin,	// minimum allowed temperature
	       Tmax;	// maximum allowed temperature

	gsl_rng * rand; // random number generator

	char * base,	// basename for animation and filename
	     * file;

	FILE * fani,	// current animation file
	     * fout;	// output file

	int * g;	// grid
} ising;

void init_ising (ising * u,
		int n, int v, int J, int mode, double h, double Tmax, double dT, double Tmin)
{
	/* we take values from the input */
	u->n      = n;
	u->h      = h;
	u->J      = J;
	u->v      = v;
	u->mode   = mode;
	u->Tmax   = Tmax;
	u->dT     = dT;
	u->Tmin   = Tmin;

	int N = n*n;
	/* prepare ground for grid initialization */
	u->g = (int *) malloc (N * sizeof (int));
	for (int i = 0; i < N-1; i++)
		u->g[i] = -1;

	/* prepare the other stuff */
	u->H  = (-1)*N*(4*J + h);
	u->m  = u->g[0]*N;

	u->base = (char *) malloc (30 * sizeof (char));
	u->file = (char *) malloc (34 * sizeof (char));

	sprintf (u->base, "J%d-N%d-h%d", u->J, u->n, (int) u->h);
	sprintf (u->file, "ISING-%s.txt", u->base);

	u->fout = fopen (u->file, "w");

	u->rand = gsl_rng_alloc (gsl_rng_mt19937);
	/* let's choose a random seed */
	FILE * frandom = fopen ("/dev/urandom", "r");

	int seed;
	fread (&seed, sizeof (int), 1, frandom);
	gsl_rng_set (u->rand, seed);

	fclose (frandom);

	/* we initialize the animation dir */
	mkdir (u->base, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	/* we will initialize T and fani later */
}

int get_spin (int * g, int N, int i, int j)
{
	return g [i*N + j];
}

void destroy_ising (ising * u)
{
	fclose (u->fout);
	free (u->base);
	free (u->file);
	free (u->g);
	free (u);
}

int step_ising (ising * u)
{
	unsigned short int i, j;

	i = (unsigned short int) u->n * gsl_rng_uniform (u->rand);
	j = (unsigned short int) u->n * gsl_rng_uniform (u->rand);

	if (i == u->n)
		i--;
	if (j == u->n)
		j--;

	int re = 0;

	int spin = get_spin (u->g, u->n, i, j),
	    s1   = get_spin (u->g, u->n, i, (j + 1)%u->n),
	    s2   = get_spin (u->g, u->n, i, (j + u->n - 1)%u->n),
	    s3   = get_spin (u->g, u->n, (i+1) % u->n, j),
	    s4   = get_spin (u->g, u->n, (i + u->n - 1) % u->n, j);

	u->dE = 2*spin * (u->J * (s1 + s2 + s3 + s4) + u->h);
	u->m += (-2)*spin;

	/* now we make "The Drawing of the Three (TM)" */
	if (u->dE < 0)
		u->g[i*u->n + j] *= (-1);

	else
	{
		float xi = gsl_rng_uniform (u->rand);

		if (xi < exp ((-1)*u->dE/(TC * u->T)))
			u->g [i*u->n + j] *= (-1);

		else
		{
			u->dE = 0;
			re    = 1;
		}
	}

	return re;
}

void update (ising * u)
{
	u->H += u->dE;
	u->E  = (1 - 1.0/u->k)*u->E  + (1.0/u->k)*u->H;
	u->E2 = (1 - 1.0/u->k)*u->E2 + (1.0/u->k)*(u->H * u->H);
	u->M  = (1 - 1.0/u->k)*u->M  + (1.0/u->k)*(1.0*u->m/(u->n * u->n));
}

void dump_regular (ising * u)
{
	double chi = (1.0/(TC * u->T)) * (1 - u->M * u->M),
	       Cv  = (1.0/pow (u->T * TC, 2)) * (u->E2 - u->E*u->E);

	fprintf (u->fout, "% 12e % 12e % 12e % 12e\n",
			u->T*TC, u->E, Cv, chi);
}

void dump_animate (ising * u)
{
	sprintf (u->file, "%s/%06d.txt", u->base, u->I);
	u->fani = fopen (u->file, "w");

	int i, j;
	for (i = 0; i <= u->n-1; i++)
	{
		for (j = 0; j <= u->n-1; j++)
			fprintf (u->fani, "% 3d", u->g[j + u->n*i]);
		fprintf (u->fani, "\n");
	}

	fclose (u->fani);
}

void dump (ising * u)
{
	dump_animate (u);
	if ( (int) (u->T * TC*500) % 10 == 0)
		dump_regular (u);
}


void reset (ising * u)
{
	u->M     = 0;
	u->E     = 0;
	u->E2    = 0;
}

// 'a' out of 'b' things are done, we started at time 'start'
void progress_bar (int a, int b, time_t * start)
{
	time_t now;
	double dsec,
	       left;
	if (start)
	{
		// we obtain current time
		time (&now);

		dsec = difftime (now, *start);
		left = (b - a) * dsec/b;

		printf ("ETA:% 3.2lfs | ", left);
	}

	int percent = 20 * a/b;
	printf ("DONE: % 6d/%d (% 3.0lf%%) [", a, b, 100.0*a/b);
	
	for (int i = 0; i <= percent-1; i++)
		printf ("=");
	printf (">");
	for (int i = percent; i <= 18; i++)
		printf (" ");
	printf ("]\n\033[F\033[J");
}

void solver (ising * u)
{
	u->T = u->Tmin;

	int max = (u->mode + 1) * (u->Tmax - u->Tmin)/u->dT;
	time_t now;

	time (&now);

	printf ("Calculating!");
	do
	{
		int r = 0;
		for (u->k = 1; u->k <= u->v; u->k++)
		{
			r += step_ising (u);
			update (u);
		}

		dump (u);
		reset (u);

		progress_bar (u->I, max, &now);
		u->I++;
		u->T += u->dT;

	} while (u->T < u->Tmax);

	if (u->mode == 1)
	{
		do
		{
			int r = 0;
			for (u->k = 1; u->k <= u->v; u->k++)
			{
				r += step_ising (u);
				update (u);
			}

			dump (u);
			reset (u);

			progress_bar (u->I, max, &now);
			u->I++;
			u->T -= u->dT;
		} while (u->T > u->Tmin);
	}
}

#endif

