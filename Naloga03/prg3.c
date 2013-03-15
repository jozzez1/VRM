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
	    arg;

	opterr = 1;
	while ((arg = getopt (argc, argv, "T:s:t:L:p:q:o:Py")) != -1)
	{
		switch (arg)
		{
			case 'T':
				T = atoi (optarg);
				break;
			case 's':
				s = atoi (optarg);
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
			default:
				abort ();
				exit (EXIT_FAILURE);
		}
	}

	if (dat == NULL)
	{
		dat = (char *) malloc (30 * sizeof (char));
		sprintf (dat, "n3-L%d-s%d-t%d",
				(int) (L*100), s, (int) (t * 100));
	}

	sprintf (file, "%s.txt", dat);
	hod * u = (hod *) malloc (sizeof (hod));
	init (u, T, L, t, q1, q2, p1, p2, file);

	printf ("Output in %s ...\n", file);
	switch (s)
	{
		case 2:
			solver2 (u);
			break;
		default:
			solver2(u);
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
