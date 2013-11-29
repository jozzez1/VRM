#include <getopt.h>
#include "hed5.h"

int main (int argc, char ** argv)
{
	// default values
	int N = 20,
	    arg;

	double dt	= 1e-2,
	       tmax	= 1000,
	       lambda	= 0.0,
	       prec	= 1e-6,
	       tau	= 1.0,
	       tdead	= 600;

	struct option longopts[] =
	{
		{ "tmax",     required_argument,        NULL,         'T' },
		{ "tdead",    required_argument,        NULL,         'd' },
		{ "prec",     required_argument,        NULL,         'p' },
		{ "dt",       required_argument,        NULL,         't' },
		{ "tau",      required_argument,        NULL,          1  },
		{ "lambda",   required_argument,        NULL,         'L' },
		{ "number",   required_argument,        NULL,         'N' },
		{ "help",     no_argument,              NULL,         'h' },
	};

	while ((arg = getopt_long (argc, argv, "T:d:p:t:L:N:h1", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'T':
				tmax = atof (optarg);
				break;
			case 'd':
				tdead = atof (optarg);
				break;
			case 'p':
				prec = atof (optarg);
				break;
			case 't':
				dt = atof (optarg);
				break;
			case 1:
				tau = atof (optarg);
				break;
			case 'L':
				lambda = atof (optarg);
				break;
			case 'N':
				N = atoi (optarg);
				break;
			case 'h':
				printf ("List of commands:\n");
				printf ("-------------------\n");
				printf ("-h / --help   ... prints this list\n");
				printf ("-T / --tmax   ... maximum time iteration\n");
				printf ("-d / --tdead  ... dead time before averaging\n");
				printf ("-t / --td     ... time step length\n");
				printf ("-L / --lambda ... lambda of the potential U(x)\n");
				printf ("-N / --number ... N = length of the chain\n");
				printf ("--tau         ... tau is the time constant for zetas\n");
				printf ("-------------------\n");
				printf ("||  defaults\n");
				printf ("-------------------\n");
				printf ("tmax	= 1000\n");
				printf ("tdead	= 600\n");
				printf ("td	= 0.01\n");
				printf ("lambda	= 0.0\n");
				printf ("N	= 20\n");
				printf ("tau	= 1\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command! Try\n%s --help\nfor list of commands.\n", argv[0]);
				exit (EXIT_FAILURE);
		}
	}
	
	hod * u = (hod *) malloc (sizeof (hod));
	init_hod (u, N, dt, tmax, lambda, prec, tau, tdead);
	solve_for_lambda (u);

	kill_hod (u);

	exit (EXIT_SUCCESS);
}

