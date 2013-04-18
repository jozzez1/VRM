// prg5.c
////////////

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed5.h"

int main (int argc, char ** argv)
{
	int N      = 10,
	    tmax   = 100,
	    tdead  = 0,
	    length = 12,
	    save   = 1,
	    arg;

	char * baseT = NULL,
	     * baseJ = NULL;
	
	double L    = 0,
	       tau  = 1,
	       TL   = 2,
	       TR   = 1,
	       prec = 1e-4,
	       h    = 1e-2;

	int dump_switch = 1;

	struct option longopts[] =
	{
		{ "TR",       required_argument,        NULL,         'R' },
		{ "TL",       required_argument,        NULL,         'L' },
		{ "tmax",     required_argument,        NULL,         'm' },
		{ "tdead",    required_argument,        NULL,         'd' },
		{ "prec",     required_argument,        NULL,         'p' },
		{ "step",     required_argument,        NULL,         'h' },
		{ "out",      required_argument,        NULL,         'o' },
		{ "tau",      required_argument,        NULL,         't' },
		{ "lambda",   required_argument,        NULL,         'l' },
		{ "number",   required_argument,        NULL,         'N' },
		{ "animate",  no_argument,              NULL,         '2' },
		{ "dump",     no_argument,              NULL,         '1' },
		{ "less",     no_argument,              NULL,         '0' },
		{ "list",     no_argument,              NULL,           1 },
	};

	while ((arg = getopt_long (argc, argv, "R:L:m:d:p:h:o:t:l:N:210", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'R':
				TR = atof (optarg);
				break;
			case 'L':
				TL = atof (optarg);
				break;
			case 'm':
				tmax = atoi (optarg);
				break;
			case 'd':
				tdead = atoi (optarg);
				break;
			case 'p':
				prec = atof (optarg);
				break;
			case 'h':
				h = atof (optarg);
				break;
			case 't':
				tau = atof (optarg);
				break;
			case 'l':
				L = atof (optarg);
				break;
			case 'N':
				N = atoi (optarg);
				break;
			case '0':
				dump_switch = 0;
				break;
			case '1':
				dump_switch = 1;
				break;
			case '2':
				dump_switch = 2;
				break;
			case 'o':
				baseJ = (char *) malloc (25 * sizeof (char));
				baseT = (char *) malloc (25 * sizeof (char));

				sprintf (baseJ, "%s-T", optarg);
				sprintf (baseT, "%s-J", optarg);

				break;
			case 1:
				printf ("List of commands:\n");
				printf ("--list,                    printf this list\n");
				printf ("--TR,       -R <1>         right-hand side temperature\n");
				printf ("--TL,       -L <2>         left-hand side temperature\n");
				printf ("--tmax,     -m <100>       maximum time iteration\n");
				printf ("--tdead,    -d <0>         dead time\n");
				printf ("--prec,     -p <1e-4>      RK4 precision\n");
				printf ("--step,     -h <1e-2>      RK4 step length\n");
				printf ("--out,      -o <format>    output file name\n");
				printf ("--tau,      -t <1>         tau value\n");
				printf ("--lambda,   -l <0>         lambda value\n");
				printf ("--number,   -N <10>        number of particles\n");
				printf ("---------------------------\n");
				printf ("--animate,  -2             animate\n");
				printf ("--dump,     -1 <default>   don't animate, just dump output\n");
				printf ("--less,     -0             only dump the end result\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s --list, for list of commands\n", argv[0]);
				abort ();
		}
	}

	/* if we have to set the name */
	if (baseJ == NULL || baseT == NULL)
	{
		baseJ = (char *) malloc (25 * sizeof (char));
		baseT = (char *) malloc (25 * sizeof (char));

		sprintf (baseJ, "J-N%d-TR%.0lf-TL%.0lf", N, TR, TL);
		sprintf (baseT, "T-N%d-TR%.0lf-TL%.0lf", N, TR, TL);
	}

	/* we initialize that sonnuvabitch */
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, N, tmax, tdead, baseT, baseJ, L, tau, TL, TR, prec, h, dump_switch);

	solver (u);
	destroy (u);

	if (dump_switch == 2)
	{
		char * command = (char *) malloc (60 * sizeof (char));
		sprintf (command, "./aniplot.sh %s %s %d %d", baseT, baseJ, length, save);

		free (baseJ);
		free (baseT);

		system (command);

		free (command);
	}

	return 0;
}

