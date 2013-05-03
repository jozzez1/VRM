// prg7.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed7.h"

int main (int argc, char ** argv)
{
	int mode       = 0,
	    jobs       = 8,
	    N          = 100,
	    M          = 100,
	    v          = 3000,
	    length     = 12,
	    arg;

	double bmax    = 5,
	       bmin    = 0,
	       db      = 0.001,
	       epsilon = 0.005,
	       lambda  = 0;

	struct option longopts[] =
	{
		{ "return",  no_argument,            &mode,        1 },
		{ "jobs",    required_argument,      NULL,       'j' },
		{ "size-M",  required_argument,      NULL,       'M' },
		{ "size-N",  required_argument,      NULL,       'N' },
		{ "wait-v",  required_argument,      NULL,       'v' },
		{ "Tmax",    required_argument,      NULL,       'T' },
		{ "Tmin",    required_argument,      NULL,       't' },
		{ "dT",      required_argument,      NULL,       'd' },
		{ "eps",     required_argument,      NULL,       'e' },
		{ "lambda",  required_argument,      NULL,       'L' },
		{ "length",  required_argument,      NULL,       'l' },
		{ "help",    no_argument,            NULL,       'h' },
		{ NULL,      0,                      NULL,         0 }
	};

	while ((arg = getopt_long (argc, argv, "j:M:N:v:T:t:d:e:L:l:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'j':
				jobs = atoi (optarg);
				break;
			case 'M':
				M    = atoi (optarg);
				break;
			case 'N':
				N    = atoi (optarg);
				break;
			case 'v':
				v    = atoi (optarg);
				break;
			case 'B':
				bmax = atof (optarg);
				break;
			case 'b':
				bmin = atof (optarg);
				break;
			case 'd':
				db   = atof (optarg);
				break;
			case 'e':
				epsilon = atof (optarg);
				break;
			case 'L':
				lambda  = atof (optarg);
				break;
			case 'l':
				length  = atoi (optarg);
				break;
			case 'h':
				printf ("===========================================\n");
				printf ("List of commands:\n");
				printf ("longopt    shortopt   default   description\n");
				printf ("===========================================\n");
				printf ("--return,               0         return back on the T scale\n");
				printf ("--jobs,    -j:        8         jobs\n");
				printf ("--size-M,  -M:        100       M\n");
				printf ("--size-N,  -N:        100       N\n");
				printf ("--wait-v,  -v:        3000      v\n");
				printf ("--Tmax,    -T:        5         max beta\n");
				printf ("--Tmin,    -t:        0         min beta\n");
				printf ("--dT,      -d:        0.001     delta beta\n");
				printf ("--eps,     -e:        0.005     epsilon\n");
				printf ("--lambda,  -L:        1         lambda\n");
				printf ("--length,  -l:        12        length in seconds\n");
				printf ("--help,    -h                   prints this list\n");
				exit (EXIT_SUCCESS);
		}
	}

	harmonic * u = (harmonic *) malloc (sizeof (harmonic));
	init_harmonic (u, jobs, N, M, v, mode, db, bmin, bmax, epsilon, lambda);

	solver (u);
	printf ("Plotting ...\n");

	int counter = 0;
	char * stuff = malloc (40 * sizeof (char));
	#pragma omp parallel shared (counter) private (stuff) num_threads (u->jobs)
	{
		#pragma omp for
		for (int k = 0; k <= u->I-1; k++)
		{
			double T = k*u->dT + u->Tmin;
			if (T > u->Tmax)
				T = 2*u->Tmax - T;

			stuff = (char *) malloc (40 * sizeof (char));
			sprintf (stuff, "./animate.sh %s %06d %d %lf",
					u->base, k, u->N, T);
			system (stuff);
			free (stuff);

			#pragma omp critical
			{
				counter++;
				progress_bar (counter, u->I-1, NULL);
			}
		}
	}

	printf ("Done!\n");

	stuff = (char *) malloc (60 * sizeof (char));
	sprintf (stuff, "./anime.sh %s %d %d %d", u->base, 1, length, u->I);
	system (stuff);
	free (stuff);

	destroy (u);

	return 0;
}

