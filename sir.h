// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <fcntl.h>
#include <limits.h>
#include <ctype.h>
#include <float.h>
#include <inttypes.h>
#include "../SFMT/SFMT.h"

#define NALLOC 0x100
#define NRND 0x100000
#define NAVG 10000000
// disease spreading macros
#define SUSCEPTIBLE 0
#define INFECTIOUS 1
#define RECOVERED 2
#define S(x) (s.state[(x)] == SUSCEPTIBLE)
#define I(x) (s.state[(x)] == INFECTIOUS)

// r n g related macro
#define RND_CHK(x) if (g.r - g.rnd >= NRND - (x)) {sfmt_fill_array32(&sfmt, g.rnd, NRND); g.r = g.rnd;}

// auxiliary macros
#define SQ(x) ((x) * (x))

typedef struct GLOBALS {
  // NETWORK SPECS
  unsigned int n, nl;
  // SIMULATION RELATED
  double r0;
  // RNG RELATED
  uint32_t *r, *rnd;
} GLOBALS;

// time steps of contacts
typedef struct NODE {
  unsigned int deg, *nb, *l;
} NODE;

// time steps of contacts
typedef struct LINK {
  unsigned int left, right;
} LINK;

// state of disease outbreak
typedef struct STATE {
  unsigned int ns, ni, nsi; // number of S, I and SI links respectively
  unsigned int *si; // a list of link indices of SI links
  unsigned int *isi; // a list of a link's index in the si list above
  unsigned int *inf; // a list of infected nodes
  unsigned short *state; // the state (S, I, R) of nodes
} STATE;

// in read.c
void read_network ();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
