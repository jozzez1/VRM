#include <hed6.h>

void update_ising (hod * u, double sigma)
{
	double w1 = exp ((-1) * (u->H + u->dE)),
	       E = ((ising *) u->phys)->E,
	       E2 = ((ising *) u->phys)->E2,
	       w = ((ising *) u->phys)->w;

	u->H += u->dE;
	((ising *) u->phys)->M -= (2.0 * sigma)/(u->n * u->n);
	((ising *) u->phys)->E  = E/(1 + w) + ((u->H * w1)/(w + w1));
	((ising *) u->phys)->E2 = E2/(1 + w) + ((u->H * u->H * w1) / (w + w1));
	((ising *) u->phys)->w += w1;
}

Bool step_ising (void * g)
{

	int i = ((hod *) g)->n * (int) gsl_rng_uniform (((hod *) g)->rand),
	    j = ((hod *) g)->n * (int) gsl_rng_uniform (((hod *) g)->rand),
	    n = ((hod *) g)->n;

	Bool success = True;

	gsl_matrix_view m = gsl_matrix_view_array (((hod *) g)->x, n, n);

	double spin = gsl_matrix_get (&m.matrix, i, j),
	       s1   = gsl_matrix_get (&m.matrix, i, (j+1) % g->n),
	       s2   = gsl_matrix_get (&m.matrix, i, (j-1+n) % n),
	       s3   = gsl_matrix_get (&m.matrix, (i + 1) % n, j),
	       s4   = gsl_matrix_get (&m.matrix, (i - 1 + n) % n, j),
	       dE   = 2*spin (h  + g->J*(s1 + s2 + s3 + s4));

	if (dE >= 0)
	{
	       xi   = gsl_rng_uniform (g->rand);

	       if (xi > exp (((-1.0) * dE)/g->T))
		       success = False;
	       else
		       gsl_matrix_set (&m.matrix, i, j, (-1)*spin);
	}

	((hod *) g)->dE = dE;

	if (success)
		update_ising ((hod *) g, spin);

	return success;
}

void dump_ising_regular (hod * u)
{
	double chi = (1.0/u->T) * (1 - pow (((ising *) u->phys)->M, 2)),
	       C   = ((ising *) u->phys)->E2;

	C -= pow (((ising *) u->phys)->E, 2);
	C *= pow(1.0/u->T, 2);

	fprintf (output->fout, "% 12e % 12e % 12e\n", u->T, chi, C);
}

void dump_ising_animate (hod * u)
{
	int i, j;

	char * filename = (char *) malloc (20 * sizeof (char));
	sprintf (filename, "%s/%05d.txt", u->basename, u->I);

	u->fani = fopen (filename, "w");

	for (i = 0; i <= u->n-1; i++)
	{
		for (j = 0; j <= u->n-1; j++)
			fprintf (u->fani, "% 2d", u->x[j + i*u->n]);
		fprintf (u->fani, "\n");
	}

	fclose (u->fani);
	free (filename);
}

void dump_ising (void * u)
{
	dump_ising_regular ((hod *) u);
	dump_ising_animate ((hod *) u);
}

void update_chain (hod * u, double s1, double s2, double s3)
{
}

Bool step_chain (void * u)
{
	int i = (int) ((hod *) g)->n * gsl_rng_uniform (((hod *) g)->rand),
	    n = ((hod *) g)->n;

	Bool success = True;

	// now we select the direction of rotation
	double costheta = (2 * gsl_rng_uniform (((hod *) g)->rand) - 1),
	       sintheta = sqrt (1 - costheta * costheta),
	       phi      = 2 * M_PI * gsl_rng_uniform (((hod *) g)->rand),
	       x1       = ((hod *) u)->x [i+1 + 0*u->n],
	       x2       = ((hod *) u)->x [i+1 + 1*u->n],
	       x3       = ((hod *) u)->x [i+1 + 2*u->n],
	       y1       = ((hod *) u)->x [i-1 + 0*u->n],
	       y2       = ((hod *) u)->x [i-1 + 1*u->n],
	       y3       = ((hod *) u)->x [i-1 + 2*u->n],
	       s1	= sin (phi) * sintheta,
	       s2	= sin (phi) * sintheta,
	       s3	= costheta,
	       J        = ((hod *) u)->J,
	       h        = ((hod *) u)->h,
	       dE       = (-1)*J*(s1*(x1 + y1) + s2*(x2 + y2) + s3*(x3 + y3)) - s3*h;

	if (dE >= 0)
	{
		xi = gsl_rng_uniform (((hod *) u)->rand);

		if (xi > exp (((-1.0) * dE)/((hod * u)->T)))
			success = False;
		else
		{
			((hod *) u)->x [i+1 + 0*u->n] += s1;
			((hod *) u)->x [i+1 + 2*u->n] += s2;
			((hod *) u)->x [i+1 + 3*u->n] += s3;
		}
	}

	if (success)
		update_chain ((hod *) u, s1, s2, s3);

	return success;
}
