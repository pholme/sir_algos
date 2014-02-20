// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Header file for the canonical algorithm.

#include <stdlib.h>
#include <math.h>
//#include <malloc.h> // needed for e.g. linux
#include "SFMT/SFMT.h"

#define NALLOC 0x100
#define NRND 0x100000
#define NAVG 100000
// disease spreading macros
#define SUSCEPTIBLE 0
#define INFECTIOUS 1
#define RECOVERED 2
#define INFECTIOUS_NEXT_TIME 3
#define S(x) (s.state[(x)] == SUSCEPTIBLE)

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
  unsigned int deg, *nb, *l;
} NODE;

// time steps of contacts
typedef struct LINK {
  unsigned int left, right;
} LINK;

// state of disease outbreak
typedef struct STATE {
  // unsigned int ns, ni, nsi; // number of S, I and SI links respectively
  // unsigned int *si; // a list of link indices of SI links
  // unsigned int *isi; // a list of a link's index in the si list above
  unsigned int *time; // time of infection
  unsigned int nchg2i, *chg2i;
  unsigned int *inf; // a list of infected nodes
  unsigned short *state; // the state (S, I, R) of nodes
} STATE;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
