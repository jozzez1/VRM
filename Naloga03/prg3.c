#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "hed3.h"

int main (int argc, char ** argv)
{
	char * dat = NULL;

	double p1 = 1,
	       p2 = 0,
	       q1 = 0,
	       q2 = 0.5,
	       L = 0,
	       t = 0.1;

	int T = 100,
	    s = 2,
	    arg;

	opterr = 1;
	while ((arg = getopt (argc, argv, "T:s:t:L:p:q:o:")) != -1)
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
				strcpy (dat, optarg);
				break;
			default:
				abort ();
		}
	}

	if (dat == NULL)
	{
		dat = (char *) malloc (20 * sizeof (char));
		sprintf (dat, "n3-L%d-s%d-t%d.txt",
				(int) (L*100), s, (int) (t * 100));
	}

	hod * u = (hod *) malloc (sizeof (hod));
	init (u, T, L, t, q1, q2, p1, p2, dat);

	free (dat);

	switch (s)
	{
		case 2:
			solver2 (u);
			break;
		default:
			solver2(u);
			break;
	}
	destroy (u);
}
