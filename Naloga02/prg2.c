// prg1.c
////////////

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "hed2.h"

int main (int argc, char ** argv)
{
	double h = 0.01,
	       L = 0,
	       a = 0,
	       t = 0.1;

	int N = 100,
	    p = 4,
	    tryrac = 1,
	    n = 0,
	    d = 0,
	    T = 100,
	    D = 8,
	    y = 0,
	    e = 0,
	    arg;

	char * dat = NULL;

	opterr = 1;

	while ((arg = getopt (argc, argv, "h:L:a:t:N:F:D:p:n:T:o:ldey")) != -1)
	{
		switch (arg)
		{
			case 'h':
				h = atof (optarg);
				break;
			case 'L':
				L = atof (optarg);
				break;
			case 'a':
				a = atof (optarg);
				break;
			case 't':
				t = atof (optarg);
				break;
			case 'N':
				N = atoi (optarg);
				break;
			case 'F':
				tryrac = atoi (optarg);
				break;
			case 'p':
				p = atoi (optarg);
				break;
			case 'n':
				n = atoi (optarg);
				break;
			case 'T':
				T = atoi (optarg);
				break;
			case 'o':
				dat = (char *) malloc (15 * sizeof (char));
				strcpy (dat, optarg);
				break;
			case 'd':
				d = 1;
				break;
			case 'y':
				y = 1;
				break;
			case 'e':
				e = 1;
				break;
			case 'D':
				D = atoi (optarg);
				break;
			case 'l':
				printf ("List of valid commands:\n\n");
				printf ("[-l] prints this list and exit\n");
				printf ("[-t] (0.1) time step\n");
				printf ("[-h] (0.01) space step\n");
				printf ("[-L] (0.0) lambda\n");
				printf ("[-a] (0.0) f(x - a)\n");
				printf ("[-N] (100) matrix rank\n");
				printf ("[-F] (1) precision test\n");
				printf ("[-T] (100) time iterations\n");
				printf ("[-p] (4) order of 2nd derivative prec.\n");
				printf ("[-n] (0) phi_n (x)\n");
				printf ("[-y] (no) save the output -- 'yes'\n");
				printf ("[-D] (8) duration of the animation in s\n");
				printf ("[-o] (format) output filename\n");
				printf ("[-e] (no) output the Energies\n");
				printf ("[-d] animation set to 'yes'\n");
				exit (EXIT_SUCCESS);
			default:
				abort ();
		}
	}

	// we initialize the program
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, h, L, a, t, N, p, d, tryrac, n, T, dat);

	// we now diagonalize the matrix
	diag_MRRR (u);
	rotate (u);

	// calculate the time stuff
	init_v (u);

	if (u->d == 1)
	{
		char * animation = (char *) malloc (40 * sizeof (char));
		sprintf (animation, "./anime.sh %s %d %d", u->dat, y, D);
		create_frames (u);
		system (animation);
	}

	else if (u->d == 0)
	{
		char * savenplot = (char *) malloc (40 * sizeof (char));
		sprintf (savenplot, "./plot.sh %s %d", u->dat, y);
		one_big_txt (u);
		system (savenplot);
	}

	if (e == 1)
	{
		char * eigen = (char *) malloc (12 * sizeof (char));
		sprintf (eigen, "Energies-N%d.txt", u->N);

		FILE * fout = fopen (eigen, "w");
		free (eigen);

		int i;
		for (i = 0; i <= u->N - 1; i++)
			fprintf (fout, "% 4d % 15.8lf\n", i, u->E[i]);

		fclose (fout);
	}
	destroy (u);

	return 0;
}
