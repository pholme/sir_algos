// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Header file for the incomplete SI link list version.

#include <stdlib.h>
#include <math.h>
#include <limits.h>
//#include <malloc.h> // needed for e.g. linux
#include "SFMT/SFMT.h"

#define NALLOC 0x100
#define NRND 0x100000
#define NAVG 100000
// disease spreading macros
#define SUSCEPTIBLE UINT_MAX
#define RECOVERED UINT_MAX - 1

// Checking if we need to generate more random numbers.
#define RND_CHK(x) if (g.r - g.rnd >= NRND - (x)) {sfmt_fill_array32(&g.sfmt, g.rnd, NRND); g.r = g.rnd;}

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
  sfmt_t sfmt;
} GLOBALS;

typedef struct NODE {
  unsigned int deg; // degree
  unsigned int ninb; // number of infected neighbors
  struct NODE **nb; // neighbors
  struct LINK **l; // adjacent links
  unsigned int t; // infection time 
} NODE;

// early is infected before late, so late can e.g. be S while early is R
typedef struct LINK {
  struct NODE *late, *early;
} LINK;

// time steps of contacts
typedef struct STATE {
  // the lchk structs contains links to check. before the SI->II stage
  // they contain SI or SR links. after the SI->II stage, they contain
  // SI or II links. they are reconstructed every time step. alternately
  // lchk1 and lchk2
  struct LINK **lchk1, **lchk2;
  struct NODE **i;
} STATE;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
