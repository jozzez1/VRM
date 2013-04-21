// prg5.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed6.h"

int main (int argc, char ** argv)
{
	int n       = 20,
	    max     = 2,
	    v       = 1000,
	    save    = 1,
	    mode    = 0,
	    length  = 12,
	    arg;

	double h    = 0,
	       J    = 1,
	       dT   = 0.01;

	struct option longopts[] =
	{
		{ "max",      required_argument,     NULL,    'm' },
		{ "magnet",   required_argument,     NULL,    'h' },
		{ "step-dT",  required_argument,     NULL,    'd' },
		{ "size-n",   required_argument,     NULL,    'n' },
		{ "force",    required_argument,     NULL,    'J' },
		{ "wait",     required_argument,     NULL,    'v' },
		{ "length",   required_argument,     NULL,    'l' },
		{ "mode",     no_argument,           NULL,    'i' },
		{ "no-save",  no_argument,           NULL,      2 },
		{ "help",     no_argument,           NULL,      1 },
	};

	while ((arg = getopt_long (argc, argv, "l:m:h:d:n:J:v:i", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'm':
				max = atoi (optarg);
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
			case 'i':
				mode = 1;
				break;

			default:
				printf ("Unknown command!\nTry %s --list, for list of commands\n", argv[0]);
				abort ();
		}
	}

	/* we initialize that sonnuvabitch */
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, mode, n, J, v, max, dT, h);

	solver (u);

	char * command = (char *) malloc (60 * sizeof (char));
	sprintf (command, "./anime.sh %s %d %d %lf", u->basename, save, length, u->dT);

	destroy (u);
	system (command);
	free (command);

	return 0;
}

