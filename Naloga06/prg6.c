// prg6.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed6.h"

int main (int argc, char ** argv)
{
	int n       = 100,
	    v       = 1000,
	    save    = 1,
	    length  = 12,
	    arg;

	double h    = 0,
	       J    = 1,
	       max  = 1.5,
	       Tmin = 0.5,
	       dT   = 0.001;

	struct option longopts[] =
	{
		{ "max",      required_argument,     NULL,    'M' },
		{ "min",      required_argument,     NULL,    'm' },
		{ "magnet",   required_argument,     NULL,    'h' },
		{ "step-dT",  required_argument,     NULL,    'd' },
		{ "size-n",   required_argument,     NULL,    'n' },
		{ "ferro",    required_argument,     NULL,    'J' },
		{ "wait",     required_argument,     NULL,    'v' },
		{ "length",   required_argument,     NULL,    'l' },
		{ "no-save",  no_argument,           NULL,      2 },
		{ "help",     no_argument,           NULL,      1 },
	};

	while ((arg = getopt_long (argc, argv, "l:M:m:h:d:n:J:v:", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'm':
				Tmin = atof (optarg);
				break;
			case 'M':
				max = atof (optarg);
				break;
			case 'h':
				h = atof (optarg);
				break;
			case 'd':
				dT = atof (optarg);
				break;
			case 'n':
				n = atoi (optarg);
				break;
			case 'J':
				J = atoi (optarg);
				if (J == 0) J = -1;
				break;
			case 'v':
				v = atoi (optarg);
				break;
			case 'l':
				length = atoi (optarg);
				break;
			case 2:
				save = 0;
				break;
			case 1:
				printf ("List of commands:\n");
				printf ("--help,          printf this list\n");
				printf ("--max,      -m   maximum allowed temperature\n");
				printf ("--length,   -l   length of the animation in seconds\n");
				printf ("--magnet,   -h   magnetic field strength\n");
				printf ("--step-dT,  -d   temperature step size\n");
				printf ("--size-n,   -n   size of our \"box\"\n");
				printf ("--force,    -J   (J = 1) ferro- or (J = -1)anti-ferromagnet\n");
				printf ("--wait,     -v   iterations we spend averaging\n");
				printf ("--no-save,       delete output after program's been finished\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s --list, for list of commands\n", argv[0]);
				abort ();
		}
	}

	/* we initialize that sonnuvabitch */
	ising * u = (ising *) malloc (sizeof (ising));
	init_ising (u, n, v, J, h, max, dT, Tmin);

	solver (u);

	char * stuff;
	#pragma omp parallel private (stuff) num_threads (8)
	{
		#pragma omp for
		for (int k = 0; k <= u->I-1; k++)
		{
			stuff = (char *) malloc (60 * sizeof (char));
			sprintf (stuff, "./animate.sh %s %06d %d %lf %lf",
					u->base, k, u->n, u->dT, u->Tmin);
			system (stuff);

			free (stuff);
		}
	}

	printf ("Finished with parallel stuff!\n");

	char * command = (char *) malloc (60 * sizeof (char));
	sprintf (command, "./anime.sh %s %d %d", u->base, save, length);

	system (command);
	free (command);
	destroy_ising (u);

	return 0;
}

