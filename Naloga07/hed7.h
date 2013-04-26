// hed7.h
////////////

#ifndef __HEADER_VRM7
#define __HEADER_VRM7

#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

typedef struct
{
	int N,
	    M,
	    jobs,
	    I,
	    mode,
	    k,
	    v;

	double E,
	       beta,
	       db,
	       bmin,
	       bmax,
	       epsilon,
	       lambda,
	       * xi;

	FILE * fout,
	     * fani;
	
	char * base,
	     * file;

	gsl_rng * rand;

	double ** c;
} harmonic;

void
init_harmonic (harmonic * h,
		int jobs,
		int N,
		int M,
		int v,
		int mode,
		double db,
		double bmin,
		double bmax,
		double epsilon,
		double lambda)
{
	h->N        = N;
	h->M        = M;
	h->v        = v;
	h->mode     = mode;
	h->jobs     = jobs;
	h->epsilon  = epsilon;
	h->lambda   = lambda;
	h->db       = db;
	h->bmax     = bmax;
	h->bmin     = bmin;
	h->beta     = bmin;

	h->k = 0;
	h->I = 0;

	/* Marsenne Twister with a REAL random seed*/
	int seed;

	FILE * frandom = fopen ("/dev/urandom", "r");
	fread (&seed, sizeof (int), 1, frandom);

	h->rand = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (h->rand, seed);

	fclose (frandom);

	/* create and populate the QM harmonic oscillator */
	h->c = (double **) malloc (h->M * sizeof (double *));
	for (int i = 0; i <= h->M-1; i++)
		h->c[i] = (double *) malloc (h->N * sizeof (double));

	h->xi = (double *) malloc (h->N * sizeof (double));

	#pragma omp parallel num_threads (h->jobs)
	{
		/* we fill that thing, can be unnormalized */
		#pragma omp for
		for (int j = 0; j <= h->M-1; j++)
		{
			for (int i = 0; i <= h->N - 1; i++)
				h->c[j][i] = gsl_ran_gaussian_ziggurat (h->rand, 1.0/h->N);
		}
	}

	/* we open our files/create animation directory */
	char * regular_dump = (char *) malloc (20 * sizeof (double));
	h->base = (char *) malloc (20 * sizeof (double));
	h->file     = (char *) malloc (30 * sizeof (double));
	
	sprintf (regular_dump,  "HARMONIC-N%d.txt", h->N);
	sprintf (h->base, "animate-N%d", h->N);

	mkdir (h->base, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	h->fout = fopen (regular_dump, "w");
	free (regular_dump);
}

void
destroy (harmonic * h)
{
	fclose (h->fout);
	free (h->base);

	for (int i = 0; i <= h->M-1; i++)
		free (h->c[i]);
	free (h->c);
	free (h);
}

double factorP (double * a, double * b, harmonic * h)
{
	double a2 = 0,
	       b2 = 0,
	       ab = 0;

	for (int i = 0; i <= h->N-1; i++)
	{
		a2 += a[i]*a[i];
		b2 += b[i]*b[i];
		ab += a[i]*b[i];
	}

	double Tp = (1.0*h->M/h->beta) * (b2 + a2 - 2 * ab),
	       Vp = (h->beta/h->M) * a2 * (0.5 + h->lambda*a2);

	return exp ((-1) * (Tp + Vp));
}

// we overwrite a with x
void overwrite (double * a, double * x, int N, int jobs)
{
	for (int i = 0; i <= N-1; i++)
		a[i] = x[i];
}

// basic Metropolis stepper function
int step_harmonic (harmonic * h)
{
	// we select the state vector at random
	int j = (int) h->M * gsl_rng_uniform (h->rand),
	    accept = 1;

	if (j == h->M)
		j--;

	// we compute the difference vector:
	// (1) we allocate it, and prepare normalization S
	double S = 0;

	// (2) so now we set the variables ...
	for (int i = 0; i <= h->N-1; i++)
	{
		h->xi[i] = gsl_ran_gaussian_ziggurat (h->rand, 1.0/h->N);
		S += h->xi[i] * h->xi[i];
	}

	// we want it to have length epsilon
	S = sqrt (S)/h->epsilon;

	// (3) we normalize it accordingly and update some vectors
	for (int i = 0; i <= h->N-1; i++)
	{
		h->xi[i] /= S;
		h->xi[i] += h->c[j][i];
	}

	// (4) now we calculate the "stuff"
	double P1p = factorP (h->c[(j+h->N-1)%h->N], h->xi, h),
	       P2p = factorP (h->xi, h->c[(j+1)%h->N], h),
	       P1  = factorP (h->c[(j-1+h->N)%h->N], h->c[j], h),
	       P2  = factorP (h->c[j], h->c[(j+1)%h->N], h),
	       p   = (P1p * P2p)/(P1 * P2);

	// now we check if we'll accept the move
	if (p < 1)
	{
		double test = gsl_rng_uniform (h->rand);
		if (test > p)
			accept = 0;
	}

	if (accept)
		overwrite (h->c[j], h->xi, h->N, h->jobs);

	return accept;
}

void dump_animate (harmonic * u)
{
	sprintf (u->file, "%s/%06d.txt", u->base, u->I);
	u->fani = fopen (u->file, "w");

	int i, j;
	for (i = 0; i <= u->N-1; i++)
	{
		for (j = 0; j <= u->M-1; j++)
			fprintf (u->fani, " %.3lf", u->c[j][i]);
		fprintf (u->fani, "\n");
	}

	fclose (u->fani);
}

void dump (harmonic * u)
{
	dump_animate (u);
}

void update (harmonic * u)
{
	return;
}

void reset (harmonic * u)
{
	return;
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

void solver (harmonic * u)
{
	u->beta = u->bmin;

	int max = (u->mode + 1) * (u->bmax - u->bmin)/u->db;
	time_t now;

	time (&now);

	printf ("Calculating ...\n");
	do
	{
		int r = 0;
		for (u->k = 1; u->k <= u->v; u->k++)
		{
			r += step_harmonic (u);
			update (u);
		}

		dump (u);
		reset (u);

		printf ("rejected: %d\t", u->v - r);
		progress_bar (u->I, max, &now);
		u->I++;
		u->beta += u->db;

	} while (u->beta < u->bmax);

	if (u->mode == 1)
	{
		do
		{
			int r = 0;
			for (u->k = 1; u->k <= u->v; u->k++)
			{
				r += step_harmonic (u);
				update (u);
			}

			dump (u);
			reset (u);

			progress_bar (u->I, max, &now);
			u->I++;
			u->beta -= u->db;
		} while (u->beta > u->bmin);
	}
}

#endif
