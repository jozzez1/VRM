// prg4.c
////////////

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "hed4.h"

int main (int argc, char ** argv)
{
	char * dat = NULL,
	     * file = (char *) malloc (30 * sizeof(char));
	
	double h = 0.01;

	int N = 3,
	    T = 0,
	    M = 100,
	    s = 4,
	    G = 10,
	    E = 0,
	    J = 0,
	    y = 0,
	    c = 0,
	    arg;

	opterr = 1;

	while ((arg = getopt (argc, argv, "o:h:N:M:s:G:tEJlc")) != -1)
	{
		switch (arg)
		{
			case 'o':
				dat = (char *) realloc ((char *) optarg, 26*sizeof(char));
				break;
			case 'h':
				h = atof (optarg);
				break;
			case 'N':
				N = atoi (optarg);
				break;
			case 't':
				T = 2;
				break;
			case 'E':
				E = 1;
				break;
			case 'J':
				J = 1;
				break;
			case 'M':
				M = atoi (optarg);
				break;
			case 's':
				s = atoi (optarg);
				break;
			case 'G':
				G = atoi (optarg);
				break;
			case 'c':
				c = 1;
				break;
			case 'y':
				y = 1;
				break;
			case 'l':
				printf ("List of commands:\n==========================\n");
				printf ("-l -- prints this list\n");
				printf ("-c -- print the connection matrix\n");
				printf ("-o -- name of the output file\n");
				printf ("-y -- save output file\n");
				printf ("-E -- calculate the average energy\n");
				printf ("-J -- calculate the the spin current\n");
				printf ("-h (0.01) set the integrator step\n");
				printf ("-N (3) set the number of qubits\n");
				printf ("-t (no) time flag\n");
				printf ("-M (100) max time iteration\n");
				printf ("-s (4) integrator precision (type)\n");
				printf ("-G (10) number of vector for averaging\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s -l for list of commands.\n", argv[0]);
				abort ();
		}
	}

	if (dat == NULL)
	{
		dat = (char *) malloc (30 * sizeof (char));

		if (T == 0)
			sprintf (dat, "zoft-G%d-N%d-E%d", G, N, E);

		else
			sprintf (dat, "corr-G%d-N%d-E%d", G, N, E);
	}

	sprintf (file, "%s.txt", dat);

	hod * u = (hod *) malloc (sizeof (hod));
	init (u, N, T, E, J, M, s, G, c, h, file);

	printf ("\nCalculating ...\n");
	simple_propagate (u);
	destroy (u);

	printf ("Done!\n");
	printf ("Output written in %s.\n", file);
	char * command = (char *) malloc (35 * sizeof (char));
	sprintf (command, "./plot4.sh %s %d %d %d", dat, y, E + T, N);
	system (command);

	free (command);
	free (file);
	free (dat);
	  
	return 0;
}

