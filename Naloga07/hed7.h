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
	       T,
	       dT,
	       Tmin,
	       Tmax,
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

void energy_start (harmonic * u)
{
	for (int i = 0; i <= u->M-1; i++)
	{
		double q1  = 0,
		       q2  = 0,
		       q12 = 0;

		for (int j = 0; j <= u->N-1; j++)
		{
			q1  += u->c[i][j] * u->c[i][j];
			q2  += u->c[(i+1)%u->N][j] * u->c[(i+1)%u->N][j];
			q12 += u->c[(i+1)%u->N][j] * u->c[i][j];
		}

		double En = (1.0/u->M) * q1*(0.5 + u->lambda * q1);
		En -= 0.5*u->T*u->T*u->M*(q1 + q2 - 2*q12);
		En /= u->N;
		En += 0.5*u->M*u->T;

		u->E = (1 - 1.0/u->k)*u->E + (1.0/u->k)*En;
	}
}

void
init_harmonic (harmonic * h,
		int jobs,
		int N,
		int M,
		int v,
		int mode,
		double dT,
		double Tmin,
		double Tmax,
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
	h->dT       = dT;
	h->Tmax     = Tmax;
	h->Tmin     = Tmin;
	h->T        = Tmin;

	h->k = 1;
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

	/* we fill that thing, can be unnormalized */
	for (int j = 0; j <= h->M-1; j++)
	{
		for (int i = 0; i <= h->N - 1; i++)
			h->c[j][i] = gsl_ran_gaussian_ziggurat (h->rand, 1.0/h->N);
	}

	/* we open our files/create animation directory */
	char * regular_dump = (char *) malloc (20 * sizeof (double));
	h->base = (char *) malloc (20 * sizeof (double));
	h->file = (char *) malloc (30 * sizeof (double));
	
	sprintf (regular_dump,  "HARMONIC-N%dL%d.txt", h->N, (int) h->lambda);
	sprintf (h->base, "video-N%d-L%d", h->N, (int) h->lambda);

	mkdir (h->base, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	h->fout = fopen (regular_dump, "w");
	free (regular_dump);

	/* we initialize the energy */
	h->E = 0;
	energy_start (h);
}

void
destroy (harmonic * h)
{
	fclose (h->fout);
	free (h->base);
	free (h->file);
	free (h->xi);

	gsl_rng_free (h->rand);

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

	double Tp = (h->M*h->T) * (b2 + a2 - 2 * ab),
	       Vp = (1.0/(h->M*h->T)) * a2 * (0.5 + h->lambda*a2);

	return exp ((-1) * (Tp + Vp));
}

// we overwrite a with x
void overwrite (double * a, double * x, int N, int jobs)
{
	for (int i = 0; i <= N-1; i++)
		a[i] = x[i];
}

void update (harmonic * u, int k)
{
	double DV  = 0,
	       DT  = 0,
	       qk  = 0,
	       q1  = 0,
	       q3  = 0,
	       qk1 = 0,
	       qk3 = 0,
	       q21 = 0,
	       q23 = 0,
	       q2  = 0;

	for (int i = 0; i <= u->N-1; i++)
	{
		q2  += u->c[k][i] * u->c[k][i];
		q3  += u->c[(k+1)%u->N][i] * u->c[(k+1)%u->N][i];
		q1  += u->c[(k+u->N-1)%u->N][i] * u->c[(k+u->N-1)%u->N][i];
		q21 += u->c[k][i] * u->c[(k+u->N-1)%u->N][i];
		q23 += u->c[k][i] * u->c[(k+1)%u->N][i];
		qk  += u->xi[i] * u->xi[i];
		qk1 += u->xi[i] * u->c[(k+u->N-1)%u->N][i];
		qk3 += u->xi[i] * u->c[(k+1)%u->N][i];
	}

	DV += qk * (0.5 + u->lambda * qk);
	DV -= q2 * (0.5 + u->lambda * q2);
	DV *= 1.0/(u->M * u->N);

	DT += qk + q3 - 2*qk3 + qk + q1 - 2*qk1;
	DT -= q2 + q3 - 2*q23 + q2 + q1 - 2*q21;
	DT *= 0.5*u->M*u->T*u->T/u->N;

	u->E += DV - DT;
}

void reset (harmonic * u)
{
	u->E = 0;
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
	{
//		update (h, j);      // broken
		overwrite (h->c[j], h->xi, h->N, h->jobs);
	}

	return accept;
}

void dump_regular (harmonic * u)
{
	fprintf (u->fout, "% 12e % 12e\n",
			u->T, u->E);
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
	dump_regular (u);
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
	printf ("DONE: % 6d/%d (% 2d%%) [", a, b, (int) round (100.0*a/b));
	
	for (int i = 0; i <= percent-1; i++)
		printf ("=");
	if (round (100.0 * a/b) != 100)
		printf (">");
	for (int i = percent; i <= 18; i++)
		printf (" ");
	printf ("]\n\033[F\033[J");
}

void solver (harmonic * u)
{
	u->T = u->Tmin;
	int max = (u->mode + 1) * (u->Tmax - u->Tmin)/u->dT;

	printf ("Calculating ...\n");
	do
	{
		int r = 0;
		for (u->k = 1; u->k <= u->v; u->k++)
		{
			r += step_harmonic (u);
			energy_start (u);
		}

		dump (u);
		reset (u);

		printf ("r: % 6d/%d  ", u->v - r, u->v);
		progress_bar (u->I-1, max, NULL);
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
				r += step_harmonic (u);
				energy_start (u);
			}

			dump (u);
			reset (u);

			printf ("r: % 6d/%d  ", u->v - r, u->v);
			progress_bar (u->I-2, max, NULL);
			u->I++;
			u->T -= u->dT;
		} while (u->T > u->Tmin);
	}
}

#endif
