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
	    C = 0,
	    y = 0,
	    c = 0,
	    arg;

	opterr = 1;

	while ((arg = getopt (argc, argv, "o:h:N:M:s:G:CEJytlc")) != -1)
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
			case 'C':
				C = 1;
				T = 1;
				break;
			case 'E':
				E = 1;
				break;
			case 'J':
				J = 1;
				T = 1;
				break;
			case 't':
				T = 1;
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
				printf ("-C -- calculate the 1st spin\n");
				printf ("-J -- calculate the the spin current\n");
				printf ("-t -- test the vector norms");
				printf ("-h (0.01) set the integrator step\n");
				printf ("-N (3) set the number of qubits\n");
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

		if (E == 1)
			sprintf (dat, "E-G%d-N%d", G, N);
		else if (C == 1)
			sprintf (dat, "C-G%d-N%d", G, N);
		else if (J == 1)
			sprintf (dat, "J-G%d-N%d", G, N);
		else
			sprintf (dat, "F-G%d-N%d", G, N);
	}

	sprintf (file, "%s.txt", dat);

	hod * u = (hod *) malloc (sizeof (hod));
	init (u, N, T, E, C, J, M, s, G, c, h, file);

	printf ("\nCalculating ...\n");
	simple_propagate (u);
	destroy (u);

	printf ("Done!\n");
	printf ("Output written in %s.\n", file);

	int mode = 0;
	if (E) mode = 1;
	else if (C) mode = 2;
	else if (J) mode = 3;

	char * command = (char *) malloc (35 * sizeof (char));
	sprintf (command, "./plot4.sh %s %d %d %d", dat, y, mode, N);
	system (command);

	free (command);
	free (file);
	free (dat);
	  
	return 0;
}

