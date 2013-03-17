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
	       E,    // Energy of the current time iteration
	       P1,   // <p_1^2>
	       P2;   // <p_2^2>

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

void updateP (hod * u)
{
	double t = u->n*u->t;
	u->P1 = (u->P1*(t - u->t) + u->x1[2]*u->x1[2]*u->t)/t;
	u->P2 = (u->P2*(t - u->t) + u->x1[3]*u->x1[3]*u->t)/t;
}

void init (hod * u, int T, double L, double t,
		double q1, double q2, double p1, double p2, char * dat)
{
	u->T = T;
	u->L = L;
	u->t = t;

	u->P1 = 0;
	u->P2 = 0;
	u->n  = 0;

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

void stepT (double c, hod * u)
{
	double k = c*u->t;

	u->x2[0] = u->x1[0] + k*u->x1[2];
	u->x2[1] = u->x1[1] + k*u->x1[3];

	u->x2[2] = u->x1[2];
	u->x2[3] = u->x1[3];

	swap (u->x1, u->x2);
}

void stepV (double c, hod * u)
{
	double k = c*u->t;

	u->x2[0] = u->x1[0];
	u->x2[1] = u->x1[1];

	u->x2[2] = u->x1[2] - k*u->x1[0]*(1 + 2*u->L*u->x1[1]*u->x2[1]);
	u->x2[3] = u->x2[3] - k*u->x1[1]*(1 + 2*u->L*u->x1[0]*u->x2[0]);

	swap (u->x1, u->x2);
}

void S1 (hod * u)
{
	stepV (1.0, u);
	stepT (1.0, u);

	Hamilton (u);
}

void S2 (double c, hod * u)
{
	stepT (c*0.5, u);
	stepV (c*1.0, u);
	stepT (c*0.5, u);

	Hamilton (u);
}

void S4 (hod * u)
{
	double x0 = (-1.0)* pow (2, 1.0/3)/(2 - pow (2, 1.0/3)),
	       x1 = 1.0 / (2 - pow (2, 1.0/3));

	S2 (x1, u);
	S2 (x0, u);
	S2 (x1, u);

	Hamilton (u);
}

void dump (hod * u)
{
	fprintf (u->fout, "% 15lf % 15lf % 15lf % 15lf % 15lf % 15lf % 15lf % 15lf\n",
			u->n*u->t, u->x1[0], u->x1[1], u->x1[2], u->x1[3],
			u->E, u->P1, u->P2);
}

void stepperS1 (hod * u)
{
	dump (u);
	do
	{
		u->n++;
		S1 (u);
		updateP (u);
		dump (u);
	} while (u->n <= u->T);
}

void stepperS2 (hod * u)
{
	dump (u);
	do
	{
		u->n++;
		S2 (1.0, u);
		updateP (u);
		dump (u);
	} while (u->n <= u->T);
}

void stepperS4 (hod * u)
{
	dump (u);
	do
	{
		u->n++;
		S4 (u);
		updateP (u);
		dump (u);
	} while (u->n <= u->T);
}

void potential (hod * u)
{
	char * dat = (char *) malloc (20*sizeof (char));
	sprintf (dat, "potential-L%d.txt", (int) (100 * u->L));
	FILE * fout = fopen (dat, "w");
	free (dat);

	double x = -1.5,
	       h = 0.01;

	do
	{
		double y = -1.5;
		do
		{
			double V = 0.5*(x*x + y*y) + u->L*x*x*y*y;
			fprintf (fout, "% 15lf % 15lf % 15lf\n", x, y, V);

			y += h;
		} while (y <= 1.5);
		fprintf (fout, "\n");

		x += h;
	} while (x <= 1.5);

	fclose (fout);
}

void destroy (hod * u)
{
	free (u->x1);
	free (u->x2);
	fclose (u->fout);
	free (u);
}
