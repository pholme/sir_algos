// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Header file for the no-state-label version.

#include <stdlib.h>
#include <math.h>
//#include <malloc.h> // needed for e.g. linux
#include "SFMT/SFMT.h"

#define NALLOC 0x100
#define NRND 0x100000
#define NAVG 100000

// Checking if we need to generate more random numbers.
#define RND_CHK(x) if (g.r - g.rnd >= NRND - (x)) {sfmt_fill_array32(&sfmt, g.rnd, NRND); g.r = g.rnd;}

// auxiliary macros
#define SQ(x) ((x) * (x))

typedef struct GLOBALS {
  // NETWORK SPECS
  unsigned int n, nl;
  // SIMULATION RELATED
  double lambda;
  unsigned int ulambda, delta;
  // RNG RELATED
  uint32_t *r, *rnd;
} GLOBALS;

typedef struct NODE {
  unsigned int deg, *nb; // degree and neighbors
  unsigned int *i; // my index in neighbor's array
  unsigned int ns; // number of susceptibles, coming first in *nb
  unsigned int t; // infection time 
} NODE;

// in read.c
void read_network ();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
