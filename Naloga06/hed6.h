// hed6.h
////////////

#ifndef __HEADER_VRM6
#define __HEADER_VRM6

#include <gsl/gsl_rand.h>
#include <gsl/gsl_matrix.h>

// our grid struct
typedef struct
{
	int * x,			// spin of the grid elements
	    n;				// grid dimensions, (n) in 1D and (n*n) in 2D case
} grid;

// observables in 2D Ising's model
typedef struct
{
	double M,			// <M>, magnetisation,
	       M2,			// <M^2>
	       E,			// <E>, energy
	       E2,			// <E^2>
	       m,			// current magentization -- for averaging
} ising;

// physical observables in Heisenberg's case
typedef struct
{
	int * C,			// spin correlations: C(r) = <s_0 * s_r>
	    n;				// Heisenberg chain length
} chain;

typedef struct
{
	grid g;				// our grid

	void * phys;			// physical observables

	double dE,			// difference in energy
	       H,			// E(sigma)
	       J,			// force parameter
	       T,			// temperature ... aww yeah!
	       h;			// magnetic field strength

	void (* update) (void *);	// function to update our variables
	void (* dump) (void);		// function to dump the data
	void (* step) (grid g);		// metropolis stepper function

	double (* test) (void *);	// function to test dE and test it against Metropolis
} hod;





#endif

