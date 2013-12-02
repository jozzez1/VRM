// hed5.h
////////////

#ifndef __HEADER_VRM5
#define __HEADER_VRM5

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

typedef struct
{
	int N;		// chain length

	double TR,	// right side temperature
	       TL,	// left side temperature
	       lambda,	// anharmonic factor lambda
	       prec,	// step precision
	       tmax,	// maximum time value
	       tdead,	// dead time -- during this we don't average
	       tau,	// time constant for zetas
	       t,	// current time
	       dt,	// maximum time step length
	       H;	// Hamiltonian energy

	double * c,	// chain itself
	       * T,	// temperature profile
	       * T2,	// profile of square temperature
	       * J,	// current
	       * J2;	// error of the current sigma = sqrt(<J^2> - <J>^2)
	
} hod;

// chain components:
///////////////////////
// 0 ... N-1  = positions of the chain
// N ... 2N-1 = chain link momenta
// 2N         = left zeta
// 2N+1       = right zeta
//

// temperature:
//////////////////////
// dT > 1:
// 	TL = 0
// 	TR = dT
// dT = 0:
// 	TL = 1
// 	TR = 1
//

double rnd (void)
{
	double re = (1.0 * (double) rand ())/RAND_MAX;
	return re;
}

// initializator
void init_hod (hod * u,
		int N,
		double dt,
		double tmax,
		double lambda,
		double prec,
		double tau,
		double tdead)
{
	u->N	= N;
	u->dt	= dt;
	u->tmax	= tmax;
	u->prec	= prec;
	u->tdead= tdead;
	u->tau	= tau;

	u->lambda= lambda;

	// temperatures on the side are hard-coded
	u->TL = 1; // values are as per instructions
	u->TR = 2; // yup ... as per instructions

	int M = 2*N + 2,
	    i;
	
	srand (500);	// we initiate a rand () with a seed value of 500
	u->c = (double *) malloc (M * sizeof (double));
	for (i = 0; i <= M-1; i++)
		u->c [i] = rnd () * 4 - 2; // random garbage

	u->T = (double *) malloc (N * sizeof (double));
	for (i = 0; i <= N-1; i++)
		u->T [i] = 0;

	u->T2= (double *) malloc (N * sizeof (double));
	for (i = 0; i <= N-1; i++)
		u->T2[i] = 0;

	u->J = (double *) malloc (N * sizeof (double));
	for (i = 0; i <= N-1; i++)
		u->J [i] = 0;

	u->J2= (double *) malloc (N * sizeof (double));
	for (i = 0; i <= N-1; i++)
		u->J2[i] = 0;
}

// we have to free the pointers at the end ...
void kill_hod (hod * u)
{
	if (u->c) free (u->c);
	if (u->T) free (u->T);
	if (u->T2)free (u->T2);
	if (u->J) free (u->J);
	if (u->J2)free (u->J2);
	if (u) free (u);
}

// the thing is f(y) = \dot{y} ... so
// f[i] = \dot{y}[i]
// y[i] = ... our chain ... ;)
int deriv (double t, const double * y, double * f, void * params)
{
	int N = ((hod *) params)->N;

	double TR	= ((hod *) params)->TR,
	       TL	= ((hod *) params)->TL,
	       lambda	= ((hod *) params)->lambda,
	       tau	= ((hod *) params)->tau;
	
	// time derivatives of positions
	int i;
	for (i = 0; i <= N-1; i++)
		f[i] = y[i + N];

	// time derivatives of momenta on the boundaries
	f[N] 	= (-1)*(2*y[0] - y[1] + lambda * pow(y[0], 3)) - y[2*N] * y[N];             // left
	f[2*N-1]= (-1)*(2*y[N-1] - y[N-2] + lambda * pow(y[N-1],3)) - y[2*N+1] * y[2*N-1]; // right

	// time derivatives of the rest of the momenta
	for (i = 1; i <= N-2; i++)
		f[i+N] = (-1)*(3 * y[i] - y[i-1] - y[i+1] + lambda * pow (y[i], 3));

	// time derivatives of the zetas
	f[2*N]		= (1.0/tau) * (pow (y[N],2) - TL);	// left
	f[2*N + 1]	= (1.0/tau) * (pow (y[2*N-1],2) - TR);	// right

	return GSL_SUCCESS;
}

// the hamiltonian ...
double energy (hod * u)
{
	double T = 0,
	       V = 0,
	       U = 0,
	       H = 0;

	int N = u->N,
	    i;

	for (i = 0; i <= N-2; i++)
	{
		T += 0.5 * pow(u->c[i+N],2);
		U += 0.5 * pow(u->c[i], 2);
		V += 0.5 * pow(u->c[i+1] - u->c[i], 2);
	}

	T += 0.5 * pow (u->c[2*N-1], 2);
	U += 0.5 * pow (u->c[N-1], 2);

	H = T + U + V;

	return H;
}

// update the average values of stuff ...
void avg (hod * u)
{
	int N = u->N,
	    i;

	for (i = 1; i <= u->N-2; i++)
	{
		u->T [i] += u->dt * pow (u->c[i+N],2);
		u->T2[i] += u->dt * pow (u->c[i+N],4);
		u->J [i] += (-0.5)*u->dt * u->c[i+N] * (u->c[i+1] - u->c[i-1]);
		u->J2[i] += u->dt * pow((-0.5) * u->c[i+N] * (u->c[i+1] - u->c[i-1]), 2);
	}

	u->T[0] += u->dt * pow (u->c[N], 2);
	u->T[N-1] += u->dt * pow (u->c[2*N-1], 2);

	u->T2[0] += u->dt * pow (u->c[N], 4);
	u->T2[N-1] += u->dt * pow (u->c[2*N-1], 4);
}

void dump (int change, hod * u)
{
	int i;
	char * dat = (char *) malloc (50 * sizeof (char));

	FILE * fout;
	switch (change)
	{
		case 0:
			// dump temperature profile
			sprintf (dat, "T-N%d-L%d.txt",
					u->N, (int) u->lambda * 100);

			fout = fopen (dat, "w");

			for (i = 0; i <= u->N-1; i++)
				fprintf (fout, "% d  % lf  % lf\n",
						i, u->T[i], sqrt(fabs(u->T2[i] - pow(u->T[i],2))));
			fprintf (fout, "\n");
			break;
		case 1:
			// dumping charge current profile
			sprintf (dat, "J-N%d-L%d.txt",
					u->N, (int) u->lambda * 100);

			fout = fopen (dat, "w");

			for (i = 0; i <= u->N-1; i++)
				fprintf (fout, "% d  % lf  % lf\n",
						i, u->J[i], sqrt(fabs(u->J2[i] - u->J[i]*u->J[i])));
			fprintf (fout, "\n");
			break;
		case 2:
			// dumping crappy energy
			sprintf (dat, "H-N%d-L%d.txt",
					u->N, (int) u->lambda * 100);

			if (u->t == 0)
				fout = fopen (dat, "w");
			else
				fout = fopen (dat, "a");

			fprintf (fout, "% lf  % e\n", u->t, u->H);
			break;
	}

	if (dat) free (dat);
	if (fout) fclose (fout);
}

void solve_for_lambda (hod * u)
{
	gsl_odeiv2_system sys = {deriv, NULL, 2*u->N + 2, u};
	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				1e-6, u->prec, 0.0);

	gsl_odeiv2_driver_set_hmax (d, u->dt);

	int N = u->tdead / u->dt,
	    i;

	double t = 0.0,
	       H = 0;

	// first we have the dead time
	for (i = 1; i <= N; i++)
	{
		double ti = i * u->dt;
		int status = gsl_odeiv2_driver_apply (d, &t, ti, u->c);

		if (status != GSL_SUCCESS)
		{
			printf ("error, return value = %d\n", status);
			break;
		}

		printf ("dead time at %lf\n", ti);

		H	= energy (u);
		u->H	= energy (u);
	}

	// now that we have thermalized our state ... so fun ensues :)
	t = 0.0;
	u->t = t;
	dump (2, u);
	N = u->tmax / u->dt;
	for (i = 1; i <= N; i++)
	{
		double ti = i * u->dt;
		int status = gsl_odeiv2_driver_apply (d, &t, ti, u->c);

		if (status != GSL_SUCCESS)
		{
			printf ("error, return value = %d\n", status);
			break;
		}

		avg (u);
		u->H = energy (u); //to see how energy changes ... if it does -- probably it does
		u->t = ti;
		printf ("t = %lf\n", u->t);
		dump (2,u);
	}

	// now we have to weight the averages with the length of the time
	// interval
	for (i = 0; i <= u->N-1; i++)
	{
		u->T [i] /= u->tmax;
		u->T2[i] /= u->tmax;
		u->J [i] /= u->tmax;
		u->J2[i] /= u->tmax;
	}

	dump (0, u);
	dump (1, u);

	if (d) gsl_odeiv2_driver_free (d);
}

void avg_J (hod * u)
{
	int i;
	double J	= 0,
	       Je2	= 0,
	       err	= 0;

	for (i = 1; i <= u->N-2; i++)
	{
		J += u->J[i];
		Je2 += u->J2[i] - u->J[i]*u->J[i];
	}
	J /= (u->N - 2);
	err = sqrt (Je2)/sqrt(u->N - 2);

	char * dat = (char *) malloc (20 * sizeof (char));
	sprintf (dat, "avg-J-L%d.txt", (int) u->lambda);

	FILE * fout = fopen (dat, "a");
	fprintf (fout, "% d  % lf %lf\n", u->N, J, err);
	
	free (dat);
	fclose (fout);
}

// add another one for scanning over lambda
void scan_through_lambda (hod * u, double dl, double lmax)
{
	do
	{
		solve_for_lambda (u);
		//
		// here comes whatever output I will think of for different lambdas ...
		//

		// now we have to iterate lambda
		u->lambda += dl;

		// let's refill the chain again with random stuff
		int i;
		for (i = 0; i <= 2*u->N + 1; i++)
			u->c [i] = 4 * rnd() - 2;
	} while (u->lambda <= lmax);
}

#endif

