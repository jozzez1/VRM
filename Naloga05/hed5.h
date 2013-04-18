// hed5.h
////////////

#ifndef __HEADER_VRM5
#define __HEADER_VRM5

#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

typedef struct
{
	int N,		// chain length
	    n,		// vector length: 2N + 2
	    t,		// current time index
	    tdead,	// time we have to wait for the problem to relaxate a bit
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
	       * avT,	// average temperature profile
	       * avJ;	// average energy curren

	// pointer to the "dump" function -- animation or no
	void (* dump) (void *);
} hod;

/* here's function for the GSL's integrator ... I give up */
int deriv (double t, const double * y, double * f, void * param)
{
	int N = ((hod *) param)->N,
	    i;

	double L   = ((hod *) param)->L,
	       tau = ((hod *) param)->tau,
	       TR  = ((hod *) param)->TR,
	       TL  = ((hod *) param)->TL;

	/* zetas on the edge */
	f[2*N]   = (1.0/tau) * (y[N]*y[N] - TL);
	f[2*N+1] = (1.0/tau) * (y[2*N-1]*y[2*N-1] - TR);

	/* q[i] derivatives are quite trivial */
	for (i = 0; i <= N-1; i++)
		f[i] = y[i+N];

	/* p[0],p[N-1] derivatives ... edge */
	f[N]     = (-1)*y[1] - y[2*N]*y[N] - y[0]*(1 + 4*L*y[0]*y[0]);
	f[2*N-1] = y[N-2] - y[2*N+1]*y[2*N-1] - y[2*N-1]*(1 + 4*L*y[2*N-1]*y[2*N-1]);

	for (i = N+1; i <= 2*N-2; i++)
		f[i] = y[i-N-1] - y[i-N+1] - y[i-N] * (1 + 4*L*y[i-N]*y[i-N]);

	return GSL_SUCCESS;
}

/* we update the average values in their respective arrays */
void update (hod * u)
{
	int i = 0,
	    t = u->t + 1,
	    N = u->N;

	u->avT[i] = ((t - 1)*1.0/t)*u->avT[i] + (1.0/t)*(0.5 * u->x[i+N]*u->x[i+N]);
	u->avJ[i] = ((t - 1)*1.0/t)*u->avJ[i] + (1.0/t) * u->x[i+1] * u->x[i+N];

	for (i = 1; i <= u->N-2; i++)
	{
		u->avT[i] = ((t - 1)*1.0/t)*u->avT[i] + (1.0/t)*(0.5 * u->x[i+N]*u->x[i+N]);
		u->avJ[i] = ((t - 1)*1.0/t)*u->avJ[i] + (1.0/t)*(u->x[i+1] - u->x[i-1]) * u->x[i+N];
	}

	u->avT[i] = ((t - 1)*1.0/t)*u->avT[i] + (1.0/t)*(0.5 * u->x[i+N]*u->x[i+N]);
	u->avJ[i] = ((t - 1)*1.0/t)*u->avJ[i] - (1.0/t) * u->x[i-1] * u->x[i+N];
}

/* function to dump these sons of bitches for animation */
void dump_animate (void * u)
{
	int i;

	char * dat1 = (char *) malloc (38 * sizeof (char)),
	     * dat2 = (char *) malloc (38 * sizeof (char));

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
		fprintf (foutJ, "% 10e", ((hod *) u)->avJ [i]); 
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
	int N = ((hod *) u)->N;

	char * dat = (char *) malloc (20 * sizeof (char));
	sprintf (dat, "TJ-final-N%d.txt", N);
	FILE * fout = fopen (dat, "w");

	for (int i = 0; i <= N-1; i++)
		fprintf (fout, "% 4d % 12e % 12e\n",
				i, ((hod *) u)->avT[i], ((hod *) u)->avJ[i]);

	fclose (fout);
	free (dat);
}

/* we finally, we solve the problem, and output the data */
void solver (hod * u)
{
	int i;

	gsl_odeiv2_system s = { &deriv, NULL, u->n, u };
	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new (&s, gsl_odeiv2_step_rk4, u->prec, u->prec, 0.0);

	double t = 0.0,
	       tnext;

	printf ("Dead time:\n");
	/* some dead steps, to reach a sufficient "weird" stat */
	for (i = 1; i <= u->tdead; i++)
	{
		printf("%d/%d\n", i, u->tdead);
		tnext = i*u->h;

		int status = gsl_odeiv2_driver_apply (d, &t, tnext, u->x);

		if (status != GSL_SUCCESS)
		{
			fprintf (stderr, "Error! Return value %d!\n", status);
			exit (1);
		}
	}

	/* now we do stuff for real */
	t = 0;
	printf ("Real time:\n");
	update (u);
	for (u->t = 1; u->t <= u->tmax; u->t++)
	{
		printf ("%d/%d\n", u->t, u->tmax);
		tnext = u->t * u->h;

		int status = gsl_odeiv2_driver_apply (d, &t, tnext, u->x);

		if (status != GSL_SUCCESS)
		{
			fprintf (stderr, "Error! Return value %d!\n", status);
			exit (1);
		}

		update (u);
		u->dump (u);
	}
	gsl_odeiv2_driver_free (d);
}

void init (hod * u,
		int N, int tmax, int tdead, char * baseT, char * baseJ,
		double L, double tau, double TL, double TR, double prec,
		double h, int dump_switch)
{
	u->N       = N;
	u->tmax    = tmax;
	u->tdead   = tdead;
	u->L       = L;
	u->tau     = tau;
	u->TL      = TL;
	u->TR      = TR;
	u->prec    = prec;
	u->h       = h;

	u->baseT = (char *) malloc (25 * sizeof (char));
	u->baseJ = (char *) malloc (25 * sizeof (char));
	strcpy (u->baseT, baseT);
	strcpy (u->baseJ, baseJ);

	if (dump_switch == 0) u->dump = &final_dump;
	if (dump_switch == 1) u->dump = &just_dump;
	if (dump_switch == 2)
	{
		u->dump = &dump_animate;

		mkdir (u->baseJ, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		mkdir (u->baseT, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}


	/* here comes the fun part ... vector initialization */
	u->n = 2 * u->N + 2;
	u->t = 0;

	u->x  = (double *) malloc (u->n * sizeof (double));
	/* we also inizialize starting points */
	u->avT = (double *) calloc (u->N, sizeof (double));
	u->avJ = (double *) calloc (u->N, sizeof (double));

	int i;
	/* the starting position is arbitrary ... let's try this one */
	/* already close to the equilibrium -- linear profile */
	for (i = 0; i <= u->N-1; i++)
		u->x [i] = 0.0;
	for (; i <= 2*u->N-1; i++)
		u->x [i] = sqrt (2 * (u->TL + (i-u->N)*(u->TR - u->TL)/(u->N-1)));
	for (; i <= u->n-1; i++)
		u->x[i] = 0.0;
}

void destroy (hod * u)
{
	free (u->baseT);
	free (u->baseJ);
	free (u->x);
	free (u->avT);
	free (u->avJ);
	free (u);
}

#endif

