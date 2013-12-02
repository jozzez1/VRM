// prg7.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed7.h"

int main (int argc, char ** argv)
{
	int N          = 100,	// chain length
	    M          = 100,	// time-cut length
	    v          = 3000,	// ... measure time?
	    arg;

	double bmax    = 20,	// maximum beta
	       bmin    = 0.01,	// minimum beta
	       db      = 0.001,	// beta step
	       epsilon = 0.05,	// MC step parameter
	       lambda  = 0;	// anharmonic potential

	struct option longopts[] =
	{
		{ "size-M",  required_argument,      NULL,       'M' },
		{ "size-N",  required_argument,      NULL,       'N' },
		{ "wait-v",  required_argument,      NULL,       'v' },
		{ "bmax",    required_argument,      NULL,       'B' },
		{ "bmin",    required_argument,      NULL,       'b' },
		{ "db",      required_argument,      NULL,       'd' },
		{ "eps",     required_argument,      NULL,       'e' },
//		{ "lambda",  required_argument,      NULL,       'L' },
		{ "lambda",  no_argument,            NULL,       'L' },
		{ "help",    no_argument,            NULL,       'h' },
		{ NULL,      0,                      NULL,         0 }
	};

	while ((arg = getopt_long (argc, argv, "M:N:v:B:b:d:e:L:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
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
				// lambda  = atof (optarg);
				lambda = 1;
				break;
//			case 'l':
//				length  = atoi (optarg);
//				break;
			case 'h':
				printf ("===========================================\n");
				printf ("List of commands:\n");
				printf ("longopt    shortopt   default   description\n");
				printf ("===========================================\n");
				printf ("--size-M,  -M:        100       M\n");
				printf ("--size-N,  -N:        100       N\n");
				printf ("--wait-v,  -v:        3000      v\n");
				printf ("--bmax,    -B:        20        max beta\n");
				printf ("--bmin,    -b:        0.01      min beta\n");
				printf ("--dT,      -d:        0.001     delta beta\n");
				printf ("--eps,     -e:        0.005     epsilon\n");
				printf ("--lambda,  -L:        1         lambda\n");
				printf ("--length,  -l:        12        length in seconds\n");
				printf ("--help,    -h                   prints this list\n");
				printf ("--no-ani,                       don't animate\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Error! Unknown parameters. Try\n");
				printf ("%s -h \tor\t %s --help\n",
						argv[0], argv[0]);
				exit (EXIT_FAILURE);
		}
	}

	hod * u = (hod *) malloc (sizeof (hod));
	init_hod (u, N, M, bmin, bmax, db, v, epsilon, lambda);

	solver (u);

	/*
	if (animate == 1)
	{
		printf ("Plotting ...\n");
		int counter = 0;
		char * stuff = malloc (60 * sizeof (char));
		for (int k = 0; k <= u->I-1; k++)
		{
			double T = k*u->dT + u->Tmin;
			if (T > u->Tmax)
				T = 2*u->Tmax - T;

			sprintf (stuff, "./animate.sh %s %06d %d %lf",
					u->base, k, u->N, T);
			system (stuff);
	
			counter++;
			progress_bar (counter, u->I-1, NULL);
		}
		printf ("Done!\n");

		sprintf (stuff, "./anime.sh %s %d %d %d", u->base, 1, length, u->I);
		system (stuff);
		free (stuff);
	}
	*/

	destroy (u);

	return 0;
}

