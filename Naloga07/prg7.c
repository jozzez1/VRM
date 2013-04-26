// prg7.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed7.h"

int main (int argc, char ** argv)
{
	int mode = 0,
	    jobs = 8,
	    N    = 100,
	    M    = 100,
	    v    = 3000;

	double bmax = 5,
	       bmin = 0,
	       db   = 0.001,
	       epsilon = 1,
	       lambda  = 0;

	harmonic * h = (harmonic *) malloc (sizeof (harmonic));
	init_harmonic (h, jobs, N, M, v, mode, db, bmin, bmax, epsilon, lambda);

	solver (h);

	destroy (h);
	return 0;
}
