// hed5.h
////////////

#ifndef __HEADER_VRM5
#define __HEADER_VRM5

// link with -lm

#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
	int N,		// chain length
	    n,		// vector length: 2N + 2
	    t,		// current time index
	    tmax;	// maximum allowed time index

	char * baseT,	// file/directory basename for temperature
	     * baseJ;	// file/directory basename for energy current
	
	double L,	// parameter lambda
	       tau,	// time parameter tau
	       TL,	// temperature of the 1st bath
	       TR,	// temperature of the 2nd bath
	       prec,	// designated step precision
	       h,	// step length
	       * x,	// x = x_t = (q_k, p_k, z_L, z_R)^T_t
	       * k1,	// k1 from RK4
	       * k2,	// k2 from RK4
	       * k3,	// k3 from RK4
	       * k4,	// k4 from RK4
	       * T,	// temporary vector for summing stuff up
	       * Y,	// teporary vector, when checking steps
	       * y,	// y = x_{t+1} = (q_k ...)^T_{t+1}
	       * avT,	// average temperature profile
	       * avJ;	// average energy curren

	// pointer to the "dump" function -- animation or no
	void (* dump) (void *);
} hod;

/* we give x, we overwrite k */
void f (double * x, double * k, hod * u)
{
	int N = u->N,
	    i;

	for (i = 0; i <= N-1; i++)
		k[i] = x[i+N];
	
	// now we do the p[0]
	k[i] = x[i-N+1] - x[i-N] * (3 + 4*u->L*x[i-N]*x[i-N]) - x[2*N]*x[i];
	++i;

	// p[i], 0 <= i <= N-1
	for (; i <= 2*N-2; i++)
		k[i] = x[i-N+1] + x[i-N-1] - x[i-N] * (3 + 4*u->L*x[i-N]*x[i-N]);

	// p[N-1]
	k[i] = x[i-N-1] - x[i-N] * (3 + 4*u->L*x[i-N]*x[i-N]) - x[2*N+1]*x[i];
	++i;

	// the baths ... z_L and z_R
	k[i] = (x[N]*x[N] - u->TL)/u->tau;
	++i;
	k[i] = (x[2*N-1]*x[2*N-1] - u->TR)/u->tau;
}

/* swap the two vectors */
void swap (double * x, double * y, hod * u)
{
	int i;
	for (i = 0; i <= u->n-1; i++)
	{
		double s = x[i];
		x[i] = y[i];
		y[i] = s;
	}
}

/* we sum two vectors with a scale h, which get's written in the temporary one */
void sum2 (double * x, double * y, double h, hod * u)
{
	int i;
	for (i = 0; i <= u->n-1; i++)
		u->T[i] = x[i] + h * y[i];
}

/* multiply vector with a scalar */
void scale (double h, double * x, hod * u)
{
	int i;
	for (i = 0; i <= u->n-1; i++)
		x[i] *= h;
}

/* norm of the difference vector */
double norm_diff (double * x, double * y, hod * u)
{
	double re = 0;
	int i;

	sum2 (x, y, (-1), u);

	for (i = 0; i <= u->n-1; i++)
		re += u->T[i] * u->T[i];

	re = sqrt (re);
	return re;
}

/* RK4 step */
void step_rk4 (double h, hod * u)
{
	f (u->x, u->k1, u);
	scale (h, u->k1, u);
	
	sum2(u->x, u->k1, 0.5, u);
	f (u->T, u->k2 , u);
	scale (h, u->k2, u);

	sum2 (u->x, u->k2, 0.5, u);
	f (u->T, u->k3, u);
	scale (h, u->k3, u);

	sum2 (u->x, u->k3, 1.0, u);
	f (u->T, u->k4, u);
	scale (h, u->k3, u);

	/* now we sum those babies up ;) */
	int i;
	for (i = 0; i <= u->n-1; i++)
		u->y[i] = u->x[i] + (u->k1[i] + u->k4[i] + 2*(u->k2[i] + u->k3[i]))/6;
}

/* we adapt the step length for precision */
void adaptive_step_rk4 (hod * u)
{
	int i = 1;

	double h = u->h/i,
	       error;

	/* we make initial step */
	step_rk4 (h, u);

	do
	{
		i++;
		h = u->h/i;

		swap (u->y, u->Y, u);

		int j;
		for (j = 0; j <= i-1; j++)
			step_rk4 (h, u);

		error = norm_diff (u->y, u->Y, u);
	} while (error >= u->prec);

	/* finally, the step was good -- we can write y into x */
	swap (u->x, u->y, u);
	u->t++;
}

/* we update the average values in their respective arrays */
void update (hod * u)
{
	int i;
	for (i = 1; i <= u->N-2; i++)
	{
		u->avT[i] *= (1 - 1.0/u->t);
		u->avT[i] += 0.5 * u->x [i+u->N] * u->x [i+u->N]/u->t;

		u->avJ[i] *= (1 - 1.0/u->t);
		u->avJ[i] += (u->x[i-1] - u->x[i+1]) * u->x[i+u->N] * 0.5/u->t;
	}

	/* we take care of those, we have previously ommited */
	/* first the temperature ... */
	u->avT[0] *= (1 - 1.0/u->t);
	u->avT[0] += 0.5 * u->x[u->N] * u->x[u->N]/u->t;

	u->avT [u->N-1] *= (1 - 1.0/u->t);
	u->avT [u->N-1] += 0.5 * u->x [2 * u->N-1] * u->x[2 * u->N-1]/u->t;

	/* ... and now the energy current */
	u->avJ[0] *= (1 - 1.0/u->t);
	u->avJ[0] += (-0.5) * u->x[1] * u->x [u->N]/u->t;

	u->avJ[u->N-1] *= (1 - 1.0/u->t);
	u->avJ[u->N-1] += 0.5 * u->x[u->N - 2] * u->x[2*u->N - 1]/u->t;
}

/* function to dump these sons of bitches for animation */
void dump_animate (void * u)
{
	int i;

	char * dat1 = (char *) malloc (30 * sizeof (char)),
	     * dat2 = (char *) malloc (30 * sizeof (char));

	sprintf (dat1, "%s/%05d.txt",
			((hod *) u)->baseT, ((hod *) u)->t);
	sprintf (dat2, "%s/%05d.txt",
			((hod *) u)->baseJ, ((hod *) u)->t);

	FILE * foutT = fopen (dat1, "w"),
	     * foutJ = fopen (dat2, "w");

	for (i = 0; i <= ((hod *) u)->N-1; i++)
	{
		fprintf (foutT, "% 5d % 10e\n", i, ((hod *) u)->avT [i]);
		fprintf (foutJ, "% 5d % 10e\n", i, ((hod *) u)->avJ [i]);
	}
	fprintf (foutT, "\n");
	fprintf (foutJ, "\n");

	fclose (foutT);
	fclose (foutJ);

	free (dat1);
	free (dat2);
}

/* in this version of dump, we just plot them in a 3d view ... we create a matrix*/
void just_dump (void * u)
{
	int i;

	char * dat1 = (char *) malloc (30 * sizeof (char)),
	     * dat2 = (char *) malloc (30 * sizeof (char));
	
	sprintf (dat1, "%s.txt", ((hod *) u)->baseT);
	sprintf (dat2, "%s.txt", ((hod *) u)->baseJ);

	FILE * foutT = fopen (dat1, "a"),
	     * foutJ = fopen (dat2, "a");
	
	for (i = 0; i <= ((hod *) u)->N-1; i++)
	{
		fprintf (foutT, "% 10e", ((hod *) u)->avT [i]);
		fprintf (foutT, "% 10e", ((hod *) u)->avJ [i]); 
	}

	fprintf (foutT, "\n");
	fprintf (foutJ, "\n");

	fclose (foutT);
	fclose (foutJ);

	free (dat1);
	free (dat2);
}

/* here we just spit out the final part */
void final_dump (void * u)
{
	int N = ((hod *) u)->N,
	    i;

	char * dat = (char *) malloc (20 * sizeof (char));
	sprintf (dat, "TJ-final-N%d.txt", N);
	FILE * fout = fopen (dat, "w");

	for (i = 0; i <= N-1; i++)
		fprintf (fout, "% 5d % 10e % 10e\n",
				((hod *) u)->t,
				((hod *) u)->avT[i],
				((hod *) u)->avJ[i]);

	fclose (fout);
	free (dat);
}

/* we finally, we solve the problem, and output the data */
void solver (hod * u)
{
	int i;
	for (i = 0; i <= u->tmax; i++)
	{
		adaptive_step_rk4 (u);
		update (u);
		u->dump (u);
	}
}

void init (hod * u,
		int N, int tmax, char * baseT, char * baseJ,
		double L, double tau, double TL, double TR,
		double prec, double h, int dump_switch)
{
	u->N    = N;
	u->tmax = tmax;
	u->L    = L;
	u->tau  = tau;
	u->TL   = TL;
	u->TR   = TR;
	u->prec = prec;
	u->h    = h;

	u->baseT = (char *) malloc (15 * sizeof (char));
	u->baseJ = (char *) malloc (15 * sizeof (char));
	strcpy (u->baseT, baseT);
	strcpy (u->baseJ, baseJ);

	if (dump_switch == 0) u->dump = &final_dump;
	if (dump_switch == 1) u->dump = &just_dump;
	if (dump_switch == 2) u->dump = &dump_animate;

	/* here comes the fun part ... vector initialization */
	u->n = 2 * u->N + 2;

	u->x   = (double *) malloc (u->n * sizeof (double));
	u->k1  = (double *) malloc (u->n * sizeof (double));
	u->k2  = (double *) malloc (u->n * sizeof (double));
	u->k3  = (double *) malloc (u->n * sizeof (double));
	u->k4  = (double *) malloc (u->n * sizeof (double));
	u->T   = (double *) malloc (u->n * sizeof (double));
	u->Y   = (double *) malloc (u->n * sizeof (double));
	u->y   = (double *) malloc (u->n * sizeof (double));

	u->avT = (double *) calloc (u->N, sizeof (double));
	u->avJ = (double *) calloc (u->N, sizeof (double));

	int i;
	/* the starting position is arbitrary ... let's try this one */
	for (i = 0; i <= u->n-1; i++)
		u->x [i] = 0.025*i*i;
}

#endif

