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
		{ "mode",    required_argument,      NULL,       'm' },
		{ "jobs",    required_argument,      NULL,       'j' },
		{ "size-M",  required_argument,      NULL,       'M' },
		{ "size-N",  required_argument,      NULL,       'N' },
		{ "wait-v",  required_argument,      NULL,       'v' },
		{ "bmax",    required_argument,      NULL,       'B' },
		{ "bmin",    required_argument,      NULL,       'b' },
		{ "db",      required_argument,      NULL,       'd' },
		{ "eps",     required_argument,      NULL,       'e' },
		{ "lambda",  required_argument,      NULL,       'L' },
		{ "length",  required_argument,      NULL,       'l' },
		{ "help",    no_argument,            NULL,       'h' },
		{ NULL,      0,                      NULL,         0 }
	};

	while ((arg = getopt_long (argc, argv, "m:j:M:N:v:B:b:d:e:L:l:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'm':
				mode = atoi (optarg);
				break;
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
				printf ("--mode,    -m:        0         mode\n");
				printf ("--jobs,    -j:        8         jobs\n");
				printf ("--size-M,  -M:        100       M\n");
				printf ("--size-N,  -N:        100       N\n");
				printf ("--wait-v,  -v:        3000      v\n");
				printf ("--bmax,    -B:        5         max beta\n");
				printf ("--bmin,    -b:        0         min beta\n");
				printf ("--db,      -d:        0.001     delta beta\n");
				printf ("--eps,     -e:        0.005     epsilon\n");
				printf ("--lambda,  -L:        1         lambda\n");
				printf ("--length,  -l:        12        length in seconds\n");
				printf ("--help,    -h                   prints this list\n");
				exit (EXIT_SUCCESS);
		}
	}

	harmonic * h = (harmonic *) malloc (sizeof (harmonic));
	init_harmonic (h, jobs, N, M, v, mode, db, bmin, bmax, epsilon, lambda);

	solver (h);
	mencoder (h, length);
	destroy (h);

	return 0;
}

