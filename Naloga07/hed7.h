// hed7.h
////////////

#ifndef __HEADER_VRM7
#define __HEADER_VRM7

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct
{
	int N,	// chain length
	    M,	// number of time frames of the chain
	    v;	//time we wait at a constant temperature, so we can average the energy

	double bmax,	// max beta
	       bmin,	// min beta
	       db,	// beta step
	       beta,	// current beta
	       eps,	// epsilon
	       L,	// quartic anharmonic parameter -- lambda
	       HpN;	// average <H>/N

	double ** q;	// chain itself
	double * x;	// random vector for stepping
} hod;

void
init_hod (hod * u,			// startup function
		int N,
		int M,
		double bmin,
		double bmax,
		double db,
		double v,
		double eps,
		double lambda)
{
	u->N 	= N;
	u->M 	= M;

	u->bmin	= bmin;
	u->bmax	= bmax;
	u->db	= db;
	u->v 	= v;

	u->eps	= eps;
	u->L	= lambda;

	// we will also initialize the random seed for our generator
	srandom (500);

	// we allocate the chain and the time frames it's in
	u->q = (double **) malloc (u->M * sizeof (double *));		// 1st index is the time-frame
	for (int j = 0; j <= u->M-1; j++)
	{
		u->q[j] = (double *) malloc (u->N * sizeof (double));	// 2nd index is the link of the chain

		// since we're in this loop, let's also initialize the chain elements
		for (int i = 0; i <= u->N-1; i++)
			u->q[j][i] = 0;
	}

	// we allocate the space for our random vector
	// and initialize its values to zero for now
	u->x = (double *) calloc (N,  sizeof (double));

	// our initial position in beta is this ...
	u->beta = u->bmin;

	// well, let's keep our intial energy a secret ... :D
	u->HpN	= 0;
}

void
destroy (hod * u)			// "destructor"
{
	for (int j = 0; j <= u->M-1; j++)
		if (u->q[j]) free (u->q[j]);
	
	if (u->q) free (u->q);
	if (u->x) free (u->x);
	
	free (u);
}

double
potential_V (hod * u, int j)
{
	// returns 1/2 (q_j^T * q_j) = 0.5 \sum_i q_ji * q_ij
	double V = 0,
	       L = 0;

	for (int i = 0; i <= u->N-1; i++)
	{
		V += pow (u->q[j % u->M][i], 2);
		L += pow (u->q[j % u->M][i], 4);	// anharmonic quartic correction
	}

	V *= 0.5;
	L *= u->L;

	V += L;

	return V;
}

double
potential_W (hod * u, int j)
{
	// not really a potential but more like "interaction" between time frames
	double W = 0;

	// periodic boundary conditions -- that's why we have "(j+1) mod M"
	for (int i = 0; i <= u->N-1; i++)
		W += pow(u->q[(j+1) % u->M][i] - u->q[j % u->M][i], 2);
	
	return W;
}

double
current_EpN (hod * u)	// E/N -- unaveraged
{
	double C   = (0.5 * u->M) / u->beta,		// constant part
	       A   = (-1) * C / (u->beta * u->N),	// weight of the W
	       B   = 1.0 / (u->M * u->N),		// weight of the V
	       S   = 0;					// sum of them all

	for (int j = 0; j <= u->M-1; j++)
		S += A * potential_W (u, j) + B * potential_V (u, j);

	S += C;
	return S;
}

double
prob_P (hod * u, int j)	// P (q_j, q_{j+1})
{
	double A = (-0.5) * (u->M / u->beta),
	       B = (-1) * (u->beta / u->M),
	       S = A * potential_W (u, j) + B * potential_V (u, j),
	       P = exp (S);

	return P;
}

void
switcheroo (int N, double * a, double * b) // we switch two vectors
{
	double c;
	for (int i = 0; i <= N-1; i++)
	{
		c	= a [i];
		a [i]	= b [i];
		b [i]	= c;
	}
}

void
replace (hod * u, int j)  // we replace the rolls of u->x and u->q[j] ... so we make a switcheroo ;)
{
	switcheroo (u->N, u->x, u->q[j]);
}

//now we need a function that does exactly one time iteration step
// how we do this? ...
//
// the general idea
// ------------------
// 1. we generate a random perturbation chain -- x
// 2. we pick a random time frame -- q[j]
// 3. calculate the P (q, x) (so that the q[j] --> q'[j] = q[j] + eps*x) (eps is a positive real number)
// 4.1 if total probability is P >= 1: we accept this move
// 4.2 if P < 1: we generate a random number in [0,1) -- y
//             a) y < P ..... we accept the move
//             b) y > P ..... we reject the move
// 
// 5. let's assume the move was accepted:
// q[j] = q'[j];
// HpN += E/N;
// we will divide the HpN with the number of time iterations at the end
//  

int
randomize_x (hod * u)		// function that creates q'[j]
{
	double N = 0;	// for normalization

	// we initialize x ... we will just use the C random number generator ... YOLO
	for (int i = 0; i <= u->N-1; i++)
	{
		u->x[i]	=  2 * (((double) random ())/RAND_MAX) - 1; // x[i] is in [-1, 1)
		N	+= u->x[i] * u->x[i];
	}

	// now we normalize this vector `x'
	for (int i = 0; i <= u->N-1; i++)
		u->x[i] = u->x[i] / sqrt (N);

	// we will make our x --> q'[j] ... so
	int j = (int) random () % u->M;
	for (int i = 0; i <= u->N-1; i++)
	{
		u->x[i] *= u->eps * sqrt(u->beta/u->bmin);
		u->x[i] += u->q[j][i];
	}

	return j;
}

int
prob_to_step (hod * u, int j)	// this function will calculate if we accept the q'[j] as the new q[j]
{
	int accept = 1;

	double Po1,	// probablity with the old q -- P(q[j-1], q[j])
	       Po2,	//          -"-      -"-     -- P(q[j], q[j+1])
	       Pn1,	// same, just replace the q[j] with the q'[j],
	       Pn2;	//                   -"-
	
	// 1st we calculate the old ones
	Po1 = prob_P (u, j);
	Po2 = prob_P (u, j+1);

	// then we switch q[j] with q'[j]
	replace (u, j);

	// now we calculate the new values
	Pn1 = prob_P (u, j);
	Pn2 = prob_P (u, j+1);

	// total probability
	double P = ((Pn1 * Pn2) / (Po1 * Po2));

	if (P < 1)
	{
		double y = ((double) random ())/RAND_MAX;
		if (y > P)
		{
			accept = 0;
			replace (u, j);	// if we reject it, we have to put the old q[j] back
		}
	}

	return accept;
}

int
one_step (hod * u)		// we only make one step with this function
{
	int j = randomize_x (u),	// we create our q'[j]
	    a = prob_to_step (u, j);	// we accept/reject it

	// we refresh the energy with either the old q[j], or the new q'[j]
	u->HpN += current_EpN (u);

	return a;
}

// now we have to move around with temperature, yes? -- yes :)
// general idea
// -----------------
// 1. make u->v steps at constant beta
// 2. on each step we sum energies together to HpN and output number of accepted moves to stdout
// 3. then we divide HpN = HpN/u->v = <H>/N
// 4. we output that number to a file as
// 	#1 beta   #2 <H>/N
// 	   ...       ...
// 5. then we reset HpN to zero and increment the beta by db
// 6. we do it, until our beta reached bmax

void
progress_bar (int a, int b)
{
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

/*
void dump_animate (hod * u, FILE * fout)
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
*/

void
solver (hod * u)
{
	// we generate a formatted filename for our output
	char * dat = (char *) malloc (60 * sizeof (char));
	sprintf (dat, "N%d-M%d-L%d.txt",
			u->N, u->M, (int) u->L);

	int part	= 0,
	    total	= (u->bmax - u->bmin) / u->db;

	// we open the file for writing, overwriting anything old
	FILE * fout = fopen (dat, "w");
	do
	{
		int a = 0;

		// we take some time to average the energy
		for (int i = 0; i <= u->v-1; i++)
			a += one_step (u);

		// we still have to divide with iterations ...
		u->HpN /= u->v;

		// we output the results
		fprintf (fout, "%lf\t%lf\t%.lf\n",
				u->beta, u->HpN, (1.0 * a)/u->v);

		//
		// or maybe even prepare them for animation -- not finished
		//

		// now we pave the way for the next step
		u->HpN = 0;
		u->beta += u->db;
		part++;

		// we note our progress
		printf ("acceptance = %.2lf\t", (1.0 * a)/u->v);
		progress_bar (part, total);
	} while (u->beta <= u->bmax);
}

#endif

