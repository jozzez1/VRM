// hed6.h
////////////

#ifndef __HEADER_VRM6
#define __HEADER_VRM6

#include <gsl/gsl_rnd.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

/* observables in 2D Ising's model */
typedef struct
{
	double M,			// <M>, magnetisation,
	       w,			// normalization of energy
	       E,			// <E>, energy
	       E2,			// <E^2>
} ising;

/* physical observables in Heisenberg's case */
typedef struct
{
	int * C,			// spin correlations: C(r) = <s_0 * s_r>
	    n;				// Heisenberg chain length
} chain;

/* the main struct ... it all depends on it */
typedef struct
{
	int * x,			// spin of the grid elements
	    n,				// grid dimensions, (n) in 1D and (n*n) in 2D case
	    v,				// annealing rate
	    I,				// current time index
	    max;			// maximum time iteration

	double h,			// magnetic field strength
	       J,			// (anti) fermomagnet parameter
	       T,			// temperature
	       dT,			// temperature step
	       dE,			// energy difference
	       H;			// E(sigma)

	FILE * fout,
	     * fani;
	
	char * basename;		// for plotting ...

	void * phys;			// physical observables -- either chain or ising

	gsl_rnd * rand;			// random number generator

	void (* step) (void *);
	void (* dump) (void *);
} hod;

/* dump function definitions -- ising */
void dump_ising_regular (hod * u);

void dump_ising_animate (hod * u);

void dump_chain_regular (hod * u);

void dump_chain_animate (hod * u);

void dump_ising (void * u);

void dump_chain (void * u);

/* stepper functions */
void update_ising (hod * u,
		double sigma);

void update_chain (hod * u,
		double sigma);

void step_ising (void * g);

void step_chain (void * g);

void solver (hod * u);

/* initializer functions */
void init_grid (int n);

void init (hod * u,
		int mode,
		int n,
		int J,
		int v,
		int max,
		double dT,
		double h);

void destroy (hod * u);

#endif

