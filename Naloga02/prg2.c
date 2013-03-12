// prg2.c
////////////

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "hed2.h"

int main (int argc, char ** argv)
{
	double h = 0,
	       L = 0,
	       a = 0,
	       t = 0.1;

	int N = 100,
	    p = 4,
	    tryrac = 1,
	    n = 0,
	    g = 3,
	    d = 0,
	    T = 100,
	    D = 8,
	    y = 0,
	    e = 0,
	    I = 0,
	    A = 0,
	    v = 0,
	    arg;

	char * dat = NULL;

	opterr = 1;

	while ((arg = getopt (argc, argv, "h:L:a:t:N:F:D:p:n:T:I:o:g:v:lAdey")) != -1)
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
			case 'g':
				g = atoi (optarg);
				break;
			case 'e':
				e = 1;
				break;
			case 'I':
				I = atoi (optarg);
				break;
			case 'A':
				A = 1;
				break;
			case 'D':
				D = atoi (optarg);
				break;
			case 'v':
				v = atoi (optarg);
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
				printf ("[-g] (0.001) threshold for good eigenvalue\n");
				printf ("[-n] (0) phi_n (x)\n");
				printf ("[-A] (no) get the best possible h for given N\n");
				printf ("[-I] (no) for calibration -- h (N)\n");
				printf ("[-v] (0) get the good eigenvalues with increment\n");
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

	if (v != 0)
	{
		if (dat == NULL)
		{
			dat = (char *) malloc (40*sizeof(char));
			sprintf (dat, "ratio-L%d.txt", (int) (100 * L));
		}

		vozi (L, a, t, g, N, p, d, tryrac, n, T, dat, v);
		printf ("Written in %s\n", dat);
		
		char * command = (char *) malloc (55 * sizeof (char));
		sprintf (command, "./plot2.sh %d %s %d", 1, dat, y);
		
		system (command);

		free (dat);
		free (command);

		exit (EXIT_SUCCESS);
	}

	// optimiziran h
	if (h == 0)
		h = 1.78850/(pow (N, 0.614827)) - N * 3.28275e-6;

	// we initialize the program
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, h, L, a, t, g, N, p, d, tryrac, n, T, dat);

	if (I == 0)
	{
		// we now diagonalize the matrix
		if (A == 0)
			diag_MRRR (u);
		else if (A == 1)
			autoh (u);
	
		if (e == 0)
		{
			rotate (u);
			init_v (u);
	
			if (u->d == 1)
			{
				char * animation = (char *) malloc (40 * sizeof (char));
				sprintf (animation, "./anime.sh %s %d %d", u->dat, y, D);
				create_frames (u);
				system (animation);
			}
		
			else
			{
				char * savenplot = (char *) malloc (40 * sizeof (char));
				sprintf (savenplot, "./plot.sh %s %d", u->dat, y);
				one_big_txt (u);
				system (savenplot);
			}
		}
	
		else
		{
			eigen_output (u);
			count_harmonic (u);
			if (y == 1)
				eigen_dump (u);
		}
	}

	else if (I != 0)
	{
		test_suite (u, I);
		printf ("Written in file h-of-N.txt\n");
		char * command = (char *) malloc (40 * sizeof (char));
		sprintf (command, "./plot2.sh 0 0 %d", y);

		system (command);
		free (command);
	}
	destroy (u);

	return 0;
}

