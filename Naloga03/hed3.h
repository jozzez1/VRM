// hed3.h
////////////

#include <math.h>
#include <stdlib.h>

typedef struct
{
	int T,       // max. time iteration
	    n;       // index of current time iteration

	double L,    // non-linear parameter Lambda
	       t,    // time step length
	       E;    // Energy of the current time iteration

	double * x1, // current time iteration
	       * x2; // next time iteration
	
	FILE * fout;
} hod;

void Hamilton (hod * u)
{
	double T = 0.5 * (u->x1[2]*u->x1[2] + u->x1[3]*u->x1[3]),
	       V = 0.5 * (u->x1[0]*u->x1[0] + u->x1[1]*u->x1[1]);
	V += u->L * u->x1[0]*u->x1[0] * u->x1[1]*u->x1[1];

	u->E = T+V;
}

void init (hod * u, int T, double L, double t,
		double q1, double q2, double p1, double p2, char * dat)
{
	u->T = T;
	u->L = L;
	u->t = t;

	u->n = 0;

	u->x1 = (double *) malloc (4 * sizeof (double));
	u->x2 = (double *) malloc (4 * sizeof (double));

	u->x1[0] = q1;
	u->x1[1] = q2;
	u->x1[2] = p1;
	u->x1[3] = p2;

	Hamilton (u);

	u->fout = fopen (dat, "w");
}

void swap (double * x1, double * x2)
{
	int i;
	for (i = 0; i <= 3; i++)
	{
		double x = x1[i];
		x1[i] = x2[i];
		x2[i] = x;
	}
}

// symmetric propagator end result
void step2 (hod * u)
{
	u->n++;

	double t2 = pow (u->t, 2);
	// q1
	u->x2[0] = (1 - t2)*u->x1[0] + 0.25*t2 + u->x1[2]*u->t;  // linear part
	u->x2[0] += t2 * u->L * u->x1[0] * u->x1[1] * u->x1[1];  // non-linear correction

	// q2
	u->x2[1] = (1 - t2)*u->x1[1] + 0.25*t2 + u->x1[3]*u->t;  // linear part
	u->x2[1] += t2 * u->L * u->x1[0] * u->x1[0] * u->x1[1];  // non-linear correction

	// p1
	u->x2[2] = (1 - t2)*u-x1[2] - 2*u->t*u->x1[0];           // linear part
	u->x2[2] -= u->L* (u->t*u->x1[0]*u->x1[1]*u->x1[1]);     // linear in t
	u->x2[2] -= u->L*t2 * (u->x1[2]*u->x1[1]*u->x1[1] + 2*u->x1[3]*u->x1[0]*u->x1[1]);

	//p2
	u->x2[3] = (1 - t2)*u-x1[3] - 2*u->t*u->x1[1];           // linear part
	u->x2[3] -= u->L* (u->t*u->x1[0]*u->x1[0]*u->x1[1]);     // linear in t
	u->x2[3] -= u->L*t2 * (u->x1[3]*u->x1[0]*u->x1[0] + 2*u->x1[2]*u->x1[0]*u->x1[1]);

	swap (u->x1, u->x2);
	Hamilton (u);
}

void dump (hod * u)
{
	fprintf (u->fout, "% 15lf % 15lf % 15lf % 15lf % 15lf % 15lf\n",
			u->n*u->t, u->x1[0], u->x1[1], u->x1[2], u->x1[3], u->E);
}

void solver2 (hod * u)
{
	do
	{
		step2 (u);
		dump (u);
	} while (u->n <= u->T);
}

void destroy (hod * u)
{
	free (u->x1);
	free (u->x2);
	fclose (u->fout);
	free (u);
}
