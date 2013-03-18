#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "hed3.h"

int main (int argc, char ** argv)
{
	char * dat = NULL,
	     * file = (char *) malloc (66 * sizeof (char)),
	     y = 'x';

	double p1 = 1,
	       p2 = 0,
	       q1 = 0,
	       q2 = 0.5,
	       L = 0,
	       t = 0.1;

	int T = 100,
	    s = 2,
	    P = 0,
	    S = 0,
	    arg;

	opterr = 1;
	while ((arg = getopt (argc, argv, "T:s:t:L:p:q:o:PylS")) != -1)
	{
		switch (arg)
		{
			case 'T':
				T = atoi (optarg);
				break;
			case 's':
				s = atoi (optarg);
				break;
			case 'S':
				S = 1;
				break;
			case 't':
				t = atof (optarg);
				break;
			case 'L':
				L = atof (optarg);
				break;
			case 'p':
				p1 = atof (optarg);
				break;
			case 'q':
				q2 = atof (optarg);
				break;
			case 'o':
				dat = (char *) realloc ((char *) optarg, 30 * sizeof (char));
				break;
			case 'P':
				P = 1;
				break;
			case 'y':
				y = 'y';
				break;
			case 'l':
				printf ("List of commands:\n");
				printf ("========================\n");
				printf ("-l <no>       -- print this list\n");
				printf ("-T (100)      -- maximum time iteration\n");
				printf ("-s (2)[1,2,4] -- stepper type\n");
				printf ("-t (0.1)      -- time step length\n");
				printf ("-L (0.0)      -- potential parameter\n");
				printf ("-p (1.0)      -- starting p1\n");
				printf ("-q (0.5)      -- starting q2\n");
				printf ("-o (format)   -- output file\n");
				printf ("-S <no>       -- scan over Lambda\n");
				printf ("-P <no>       -- calculate potential too\n");
				printf ("-y <no>       -- save output\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s -l for help\n", argv[0]);
				abort ();
				exit (EXIT_FAILURE);
		}
	}

	if (S == 1)
	{
		dat = (char *) malloc (30 * sizeof (char));
		file= (char *) malloc (34 * sizeof (char));

		sprintf (dat, "scan-t%d", (int) (t * 100));
		sprintf (file,"%s.txt", dat);

		hod * u = (hod *) malloc (sizeof (hod));
		init (u, T, 0, t, q1, q2, p1, p2, s, file);
		free (file);

		Lscan (u);
		destroy (u);
		printf ("Output in %s.txt ...\n", dat);

		int save = 0;
		if (y == 'y')
			save = 1;


		char * command = (char *) malloc (50*sizeof (char));
		sprintf (command, "./plot3L.sh %s %d",
				dat, save);

		system (command);

		free (command);
		free (dat);

		exit (EXIT_SUCCESS);
	}

	if (dat == NULL)
	{
		dat = (char *) malloc (30 * sizeof (char));
		sprintf (dat, "n3-L%d-s%d-t%d",
				(int) (L*100), s, (int) (t * 100));
	}

	sprintf (file, "%s.txt", dat);
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, T, L, t, q1, q2, p1, p2, s, file);

	printf ("Output in %s ...\n", file);
	switch (s)
	{
		case 1:
			stepperS1 (u);
			break;
		case 2:
			stepperS2 (u);
			break;
		case 4:
			stepperS4 (u);
			break;
		default:
			stepperS2(u);
			break;
	}

	if (P == 1)
		potential (u);

	sprintf (file, "./plot3.sh %s %d %c", dat, (int) (100 * L), y);
	system (file);

	destroy (u);
	free (file);
	free (dat);

	return 0;
}

