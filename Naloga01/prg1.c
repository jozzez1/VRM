// prg1.c
////////////

#include <stdio.h>
#include <string.h>
#include <ctype.h> // for integer conversion
#include <stdlib.h>
#include <sys/stat.h>  // mkdir ()
#include <sys/types.h> // mkdir ()
#include <unistd.h>
#include "hed1.h"

int main (int argc, char ** argv)
{
	double t = 0.1,
	       h = 0,
	       a = 0.0,
	       Lambda = 0.0;

	int N = 100,
	    M = 100,
	    s = 0,
	    y = 0,
	    T = 12,
	    dflag = 0,
	    pet = 0,
	    arg;

	char * dat = NULL;

	opterr = 1;

	while ((arg = getopt (argc, argv, "t:h:a:L:N:M:s:o:T:ldyp")) != -1)
	{
	       switch (arg)
	       {
		       case 't':
			       t = atof (optarg);
			       break;

		       case 'h':
			       h = atof (optarg);
			       break;

		       case 'a':
			       a = atof (optarg);
			       break;

		       case 'L':
			       Lambda = atof (optarg);
			       break;

		       case 'N':
			       N = atoi (optarg);
			       break;

		       case 'M':
			       M = atoi (optarg);
			       break;

		       case 's':
			       s = atoi (optarg);
			       break;

		       case 'y':
			       y = 1;
			       break;

		       case 'T':
			       T = atoi (optarg);
			       break;

		       case 'p':
			       pet = 1;
		               break;

		       case 'o':
			       dat = (char *) malloc (40 * sizeof (char));
			       strcpy (dat, optarg);
			       break;

		       case 'd':
			       dflag = 1;
			       break;

		       case 'l':
			       printf ("List of valid commands:\n\n");
			       printf (" -l <dbl> -- prints this list\n");
			       printf (" -t <dbl> (0.1) -- time granularity - tau\n");
			       printf (" -h <dbl> (3*t) -- space granularity - h\n");
			       printf (" -a <dbl> (0.0) -- start excentricity - alpha\n");
			       printf (" -L <dbl> (0.0) -- potential parameter - Lambda\n");
			       printf (" -N <int> (100) -- number of time-grains\n");
			       printf (" -M <int> (100) -- number of space-grains\n");
			       printf (" -s [0,1,2] (0) -- eigen function mode\n");
			       printf (" -d (no) -- animation mode (opposite is plotting)\n");
			       printf (" -y (no) -- save the end file");
			       printf (" -T <int> (12) -- animation length in seconds");
			       printf (" -p (no) -- use five-diagonal matrix instead");
			       printf (" -o <char *> -- name of output file\n\n");
			       exit (EXIT_SUCCESS);

		       default:
			       abort ();
	       }
	}

	if (h == 0)
		h = 3*t;
	// we arrange the filename
	
	if (dat == NULL)
	{
		dat = (char *) malloc (40 * sizeof (char));
		sprintf (dat, "phi%d_N%d_M%d_%d_%.0lf",
				s, N, M, (int) (100 * Lambda), a);
	}

	char * file1 = (char *) malloc (44 * sizeof(char));
	sprintf (file1, "%s.txt", dat);

	// initialize the "walker"
	hod * u = (hod *) malloc (sizeof(hod));

	// allocate internal varibales
	init (u, h, t, a, Lambda, N, M, s, pet);

	strcpy (u->dat, dat);

	free (dat);
	// solve the problem ...
	if (dflag == 0)
	{
		u->fout = fopen (file1, "w");

		char * savenplot = (char *) malloc (40 * sizeof(char));
		sprintf (savenplot, "./plot.sh %s %d", u->dat, y);

		solver (u);
		printf ("Writing in %s\n", file1);

		fclose (u->fout);
		// deallocate the space
		destroy (u);

		system (savenplot);
		free (savenplot);
	}

	if (dflag == 1)
	{
		char * animation = (char *) malloc (60 * sizeof(char));
		sprintf (animation, "./anime.sh %s %d %d", u->dat, y, T);

		int status;
		status = mkdir (u->dat, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		shitter (u);
		printf ("Created the %s directory\n", u->dat);
		// destroy the sonnuvabitch
		destroy (u);

		// call the animation zsh script
		system (animation);
		free (animation);
	}

	return 0;
}

