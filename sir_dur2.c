// The no-state-label version of the SIR model.

#include "sir_dur2.h"

struct GLOBALS g;
NODE *n;
unsigned int *inf;
sfmt_t sfmt;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void read_network (char *fname) {
    unsigned int i, me, you;
    FILE *fp;
    
    fp = fopen(fname, "r");
    if (!fp) {
        fprintf(stderr, "can't open %s\n", fname);
        exit(1);
    }
    g.n = 0;
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
        if (g.n < me) g.n = me;
        if (g.n < you) g.n = you;
    }
    rewind(fp);
    g.n++;
    
    n = (NODE *) malloc(g.n * sizeof(NODE));
    for (i = 0; i < g.n; i++) n[i].deg = 0;
    
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
        n[me].deg++;
        n[you].deg++;
    }
    
    for (i = 0; i < g.n; i++) {
        n[i].nb = (unsigned int *) malloc(n[i].deg * sizeof(unsigned int));
        n[i].i = (unsigned int *) malloc(n[i].deg * sizeof(unsigned int));
        n[i].deg = 0;
    }
    rewind(fp);
    
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
        n[me].i[n[me].deg] = n[you].deg;
        n[you].i[n[you].deg] = n[me].deg;
        n[me].nb[n[me].deg++] = you;
        n[you].nb[n[you].deg++] = me;
    }
    fclose(fp);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// me is doing S->I, whose neighbors need to adjust their arrays &
// move the infectious to behind the susceptibles

void swapdown (unsigned int me) {
    unsigned int you, other, me_iyou, you_iother, you_ime, other_iyou;
    
    for (me_iyou = 0; me_iyou < n[me].deg; me_iyou++) {
        you = n[me].nb[me_iyou];
        
        // other is the last of you's susceptibles
        you_iother = --n[you].ns;
        other = n[you].nb[you_iother];
        
        // to avoid unnecessary swapping, not needed but speeds it up
        if (other != me) {
            you_ime = n[me].i[me_iyou];
            
            // you's index among other's neighbors
            other_iyou = n[you].i[you_iother];
            
            // other's index among you's neighbors will be me's old index
            n[other].i[other_iyou] = you_ime;
            
            // me's index among you's neighbors will be other's old index
            n[me].i[me_iyou] = you_iother;
            
            // do the neighbor swapping
            n[you].nb[you_iother] = me;
            n[you].nb[you_ime] = other;
            
            // swap the corresponding indices
            n[you].i[you_iother] = me_iyou;
            n[you].i[you_ime] = other_iyou;
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double sir () {
    unsigned int i, me, ime, you, now, iold = 0, inew = 1, inew_now;
    
    // initialize; infect a source
    // get a random seed
    RND_CHK(1);
    me = *(g.r++) % g.n;
    
    // add me as infectious
    for (i = 0; i < g.n; i++) n[i].ns = n[i].deg;
    inf[0] = me;
    n[me].t = 0;
    swapdown(me);
    
    // run the updates
    for (now = 1; ; now++) {
        // I -> R  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for ( ; iold < inew; iold++)
            if (n[inf[iold]].t + g.delta >= now) break;
        
        // if the outbreak has died, break the iterations
        if (inew <= iold) return inew / (double) g.n;
        
        // SI -> II  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // go over the infected nodes me
        inew_now = inew; // so we don't spread the disease further than one step
        for (ime = iold; ime < inew_now; ime++) {
            me = inf[ime];
            // check me's susceptible neighbors
            RND_CHK(n[me].ns);
            for (i = 0; i < n[me].ns; ) {
                you = n[me].nb[i];
                if (*(g.r++) < g.ulambda) {
                    // infect you
                    inf[inew++] = you;
                    n[you].t = now;
                    swapdown(you);
                } else i++; // n[me].ns decreases in swapdown, needn't increment i
            }
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
        fprintf(stderr, "usage: sir_dur2 [rng state file] [nwk file] [lambda] [delta]\n");
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
    // for the transmission probability in units of a 32-bit int
    g.ulambda = (unsigned int) rint(UINT32_MAX * g.lambda);
    g.delta = atoi(argv[4]);
    
    // alloc
    inf = (unsigned int *) malloc(g.n * sizeof(unsigned int));
    
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
    for (i = 0; i < g.n; i++) {free(n[i].nb); free(n[i].i);}
    free(n); free(g.rnd); free(inf);
    
    return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
