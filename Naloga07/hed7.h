// hed7.h
////////////

#ifndef __HEADER_VRM7
#define __HEADER_VRM7

#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

typedef struct
{
	int N,
	    M,
	    jobs,
	    I,
	    k,
	    v;

	double E,
	       beta,
	       epsilon,
	       P,
	       lambda;

	FILE * fout,
	     * fani;
	
	char * base_dir;

	gsl_rng * rand;

	double ** c;
} harmonic;

void
init_harmonic (harmonic * h,
		int jobs,
		int N,
		int M,
		int v,
		double beta,
		double epsilon,
		double lambda)
{
	h->N        = N;
	h->M        = M;
	h->v        = v;
	h->jobs     = jobs;
	h->beta     = beta;
	h->epsilon  = epsilon;
	h->lambda   = lambda;

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
	h->base_dir = (char *) malloc (20 * sizeof (double));
	
	sprintf (regular_dump,  "HARMONIC-N%d.txt", h->N);
	sprintf (h->base_dir, "animate-N%d", h->N);

	mkdir (h->base_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	h->fout = fopen (regular_dump, "w");
	free (regular_dump);
}

void
destroy (harmonic * h)
{
	fclose (h->fout);
	free (h->base_dir);

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

	#pragma omp parallel shared (a2, b2, ab) num_threads (h->jobs)
	{
		#pragma omp for reduction (+:a2, b2, ab)
		for (int i = 0; i <= h->N-1; i++)
		{
			a2 += a[i]*a[i];
			b2 += b[i]*b[i];
			ab += a[i]*b[i];
		}
	}

	double Tp = (1.0*h->M/h->beta) * (b2 + a2 - 2 * ab),
	       Vp = (h->beta/h->M) * a2 * (0.5 + h->lambda*a2);

	return exp ((-1) * (Tp + Vp));
}

// we overwrite a with x
void overwrite (double * a, double * x, int N, int jobs)
{
	#pragma omp parallel num_threads (jobs)
	{
		#pragma omp for
		for (int i = 0; i <= N-1; i++)
			a[i] = x[i];
	}
}

// basic Metropolis stepper function
int harmonic_step (harmonic * h)
{
	// we select the state vector at random
	int j = (int) h->M * gsl_rng_uniform (h->rand);

	if (j == h->M)
		j--;

	// we compute the difference vector:
	// (1) we allocate it, and prepare normalization S
	double * xi = (double *) malloc (h->N * sizeof (double)),
	       S    = 0;

	#pragma omp parallel shared (xi, S) num_threads (h->jobs)
	{
		// (2) so now we set the variables ...
		#pragma omp for reduction (+:S)
		for (int i = 0; i <= h->N-1; i++)
		{
			xi[i] = gsl_ran_gaussian_ziggurat (h->rand, 1.0/h->N);
			S += xi[i] * xi[i];
		}

		// we want it to have length epsilon
		S = sqrt (S)/h->epsilon;

		// (3) we normalize it accordingly and update some vectors
		#pragma omp for
		for (int i = 0; i <= h->N-1; i++)
		{
			xi[i] /= S;
			xi[i] += h->c[j][i];
		}
	}

	// (4) now we calculate the "stuff"
	double P1p = factorP (h->c[j-1], xi, h),
	       P2p = factorP (xi, h->c[j+1], h),
	       P1  = factorP (h->c[j-1], h->c[j], h),
	       P2  = factorP (h->c[j], h->c[j+1], h),
	       p   = (P1p * P2p)/(P1 * P2);

	// now we check if we'll accept the move
	int accept = 1;

	if (p < 1)
	{
		double test = gsl_rng_uniform (h->rand);
		if (test > p)
			accept = 0;
	}

	if (accept)
		overwrite (h->c[j], xi, h->N, h->jobs);

	return accept;
}

// we concatenate these steps to obtain the solver
void solver (harmonic * h)
{

}

#endif
