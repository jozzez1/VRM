// hed3.h
////////////

#ifndef __HEADER_VRM3
#define __HEADER_VRM3

#include <math.h>
#include <stdlib.h>

//RK4 stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

typedef struct
{
	int T,       // max. time iteration
	    n,       // index of current time iteration
	    s;

	double L,    // non-linear parameter Lambda
	       t,    // time step length
	       E,    // Energy of the current time iteration
	       E0,   // starting energy
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
		double q1, double q2, double p1, double p2, int s, char * dat)
{
	u->T = T;
	u->L = L;
	u->t = t;
	u->s = s;

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
	u->E0 = u->E;

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
	fprintf (u->fout, "% 15lf % 15lf % 15lf % 15lf % 15lf % 15lf % 15lf % 15lf % 15.8e\n",
			u->n*u->t, u->x1[0], u->x1[1], u->x1[2], u->x1[3],
			u->E, u->P1, u->P2, u->E - u->E0);
}

void dumpL (hod * u)
{
	fprintf (u->fout, "% 15lf % 15lf % 15lf\n", u->L, u->P1, u->P2);
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

void Lstepper (hod * u, double L)
{
	u->L = L;

	u->P1 = 0;
	u->P2 = 0;
	u->n  = 0;

	u->x1[0] = 0;
	u->x1[1] = 0.5;
	u->x1[2] = 1.0;
	u->x1[3] = 0.0;

	Hamilton (u);
	do
	{
		u->n++;

		switch (u->s)
		{
			case 1:
				S1 (u);
				break;
			case 2:
				S2 (1.0, u);
				break;
			case 4:
				S4 (u);
				break;
			default:
				S4 (u);
				break;
		}
		updateP (u);
	} while (u->n <= u->T);
}

void Lscan (hod * u)
{
	double L  = 30,
	       dL = 0.01;

	printf ("Lambda = % 6.2lf\n", u->L);
	Lstepper (u, 0.0);
	dumpL (u);
	do
	{
		printf ("Lambda = % 6.2lf\n", u->L + dL);
		Lstepper (u, u->L + dL);
		dumpL (u);
	} while (u->L <= L);
}

void potential (hod * u)
{
	char * dat = (char *) malloc (20*sizeof (char));
	sprintf (dat, "potential-L%d.txt", (int) (100 * u->L));
	FILE * fout = fopen (dat, "w");
	free (dat);

	double x = -1.5,
	       h = 0.05;

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

int funcRK4 (double t, const double * y, double * f, void * params)
{
	double L = * (double *) params;
	f[0] = y[2];
	f[1] = y[3];
	f[2] = (-1)*y[0]*(1 + L*y[1]*y[1]);
	f[3] = (-1)*y[1]*(1 + L*y[0]*y[0]);
	return GSL_SUCCESS;
}

void RK4 (hod * u)
{
	gsl_odeiv2_system sys = {funcRK4, NULL, 4, &u->L};
	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				u->t, 1e-8, 1e-8);

	int s;

	double t = u->n*u->t;
	dump (u);
	do
	{
		u->n++;
		s = gsl_odeiv2_driver_apply_fixed_step (d, &t, u->t, 1000, u->x1);

		if (s != GSL_SUCCESS)
		{
			printf ("Driver error %d\n", s);
			break;
		}

		Hamilton (u);
		updateP (u);
		dump (u);
	} while (u->n <= u->T);

	gsl_odeiv2_driver_free (d);
}

void destroy (hod * u)
{
	free (u->x1);
	free (u->x2);
	fclose (u->fout);
	free (u);
}

#endif
