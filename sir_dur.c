// The canonical algorithm

#include "sir_dur.h"

GLOBALS g;
LINK *l;
NODE *n;
STATE s;
sfmt_t sfmt;
unsigned int allocl;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void add_link (unsigned int me, unsigned int you) {
    
    if (allocl <= g.nl) {
        allocl = g.nl + NALLOC;
        l = (LINK *) realloc(l, allocl * sizeof(LINK));
    }
    l[g.nl].left = me;
    l[g.nl++].right = you;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void read_network (char *fname) {
    unsigned int i, me, you;
    FILE *fp;
    
    fp = fopen(fname, "r");
    if (!fp) {
        fprintf(stderr, "can't open %s\n", fname);
        exit(1);
    }
    allocl = g.nl = 0;
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) add_link(me, you);
    fclose(fp);
    
    l = (LINK *) realloc(l, g.nl * sizeof(LINK));
    
    g.n = 0;
    for (i = 0; i < g.nl; i++) {
        if (g.n < l[i].left) g.n = l[i].left;
        if (g.n < l[i].right) g.n = l[i].right;
    }
    g.n++;
    
    n = (NODE *) malloc(g.n * sizeof(NODE));
    for (i = 0; i < g.n; i++) n[i].deg = 0;
    
    for (i = 0; i < g.nl; i++) {
        n[l[i].left].deg++;
        n[l[i].right].deg++;
    }
    
    for (i = 0; i < g.n; i++) {
        n[i].nb = (unsigned int *) malloc(n[i].deg * sizeof(unsigned int));
        n[i].l = (unsigned int *) malloc(n[i].deg * sizeof(unsigned int));
        n[i].deg = 0;
    }
    
    for (i = 0; i < g.nl; i++) {
        me = l[i].left;
        you = l[i].right;
        n[me].nb[n[me].deg] = you;
        n[you].nb[n[you].deg] = me;
        n[me].l[n[me].deg++] = i;
        n[you].l[n[you].deg++] = i;
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double sir () {
    unsigned int i, me, ime, you, now, iold = 0, inew = 1;
    
    // initialize; infect a source
    // get a random seed
    RND_CHK(1);
    me = *(g.r++) % g.n;
    // initialize the states
    for (i = 0; i < g.n; i++) s.state[i] = SUSCEPTIBLE;
    // add me as infectious
    s.state[me] = INFECTIOUS;
    s.inf[0] = me;
    s.time[me] = 0;
    
    // run the updates
    for (now = 1; ; now++) {
        // I -> R
        for ( ; iold < inew; iold++) {
            me = s.inf[iold];
            if (s.time[me] + g.delta < now) s.state[me] = RECOVERED;
            else break;
        }
        
        // if the outbreak has died, break the iterations
        if (inew <= iold) return inew / (double) g.n;
        
        // SI -> II
        for (ime = iold, s.nchg2i = 0; ime < inew; ime++) {
            me = s.inf[ime];
            // Scanning the neighborhood of the infectious node me.
            for (i = 0; i < n[me].deg; i++) {
                you = n[me].nb[i];
                // If the neighbor you is susceptible.
                if (S(you)) {
                    RND_CHK(1);
                    // With a probability lambda . .
                    if (*(g.r++) < g.ulambda) {
                        // infect you
                        s.chg2i[s.nchg2i++] = you;
                        s.state[you] = INFECTIOUS_NEXT_TIME;
                    }
                }
            }
        }
        
        // chg the ones to be infectious next time step
        for (i = 0; i < s.nchg2i; i++) {
            me = s.chg2i[i];
            s.inf[inew++] = me;
            s.state[me] = INFECTIOUS;
            s.time[me] = now;
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input

int main (int argc, char *argv[]) {
    unsigned int i;
    double dtmp, s1 = 0.0, s2 = 0.0;
    FILE *fp;
    
    if (argc != 5) {
        fprintf(stderr, "usage: sir_dur [rng state file] [nwk file] [lambda] [delta]\n");
        return 1;
    }
    
    // read state or initialize RNG
    // on linux malloc should be changed to memealign and malloc.h loaded
    // g.rnd = (uint32_t *) memalign(16, NRND * sizeof(uint32_t));
    g.rnd = (uint32_t *) malloc(NRND * sizeof(uint32_t));
    
    
    fp = fopen(argv[1], "rb");
    if (!fp) {
        fprintf(stderr, "can't open %s\n", argv[1]);
        return 1;
    }
    if (1 != fread(&sfmt, sizeof(sfmt_t), 1, fp)) {
        fprintf(stderr, "error reading rng state\n");
        return 1;
    }
    fclose(fp);
    sfmt_fill_array32(&sfmt, g.rnd, NRND);
    g.r = g.rnd;
    // read state or initialize RNG, end
    
    read_network(argv[2]);
    g.lambda = atof(argv[3]);
    g.ulambda = (unsigned int) rint(UINT32_MAX * g.lambda);
    g.delta = atoi(argv[4]);
    
    // alloc
    s.chg2i = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    s.inf = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    s.time = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    s.state = (unsigned short *) malloc(g.n * sizeof(unsigned short));
    
    // disease simulation
    for (i = 0; i < NAVG; i++) {
        dtmp = sir();
        s1 += dtmp;
        s2 += SQ(dtmp);
    }
    s1 /= NAVG; s2 /= NAVG;
    printf("%g %g\n", s1, sqrt((s2 - SQ(s1)) / (NAVG - 1.0)));
    
    // store rng state
    fp = fopen(argv[1], "wb");
    if (!fp) {
        fprintf(stderr, "can't open %s\n", argv[1]);
        return 1;
    }
    if (1 != fwrite(&sfmt, sizeof(sfmt_t), 1, fp)) {
        fprintf(stderr, "error writing rng state\n");
        return 1;
    }
    fclose(fp);
    
    // cleaning
    for (i = 0; i < g.n; i++) {
        free(n[i].nb);
        free(n[i].l);
    }
    free(n); free(l); free(g.rnd);
    free(s.chg2i); free(s.state); free(s.inf); free(s.time);
    
    return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
