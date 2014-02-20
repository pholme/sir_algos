// The incomplete SI link list version.
#include "sir_dur3.h"

struct GLOBALS g;
struct LINK *l;
struct NODE *n;
struct STATE s;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Helper routine for the network reading routine.

void add_link (unsigned int me, unsigned int you) {
    
    n[me].l[n[me].deg] = n[you].l[n[you].deg] = l + g.nl;
    l[g.nl].early = n[me].nb[n[me].deg++] = n + you;
    l[g.nl++].late = n[you].nb[n[you].deg++] = n + me;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Reading the network, it should be a white space separated edge list with
// node indices from 0 to N - 1

void read_network (char *fname) {
    unsigned int i, me, you;
    FILE *fp;
    
    fp = fopen(fname, "r");
    if (!fp) {
        fprintf(stderr, "can't open %s\n", fname);
        exit(1);
    }
    g.n = g.nl = 0;
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
        if (g.n < me) g.n = me;
        if (g.n < you) g.n = you;
        g.nl++;
    }
    rewind(fp);
    g.n++;
    
    n = (NODE *) malloc(g.n * sizeof(NODE));
    l = (LINK *) malloc(g.nl * sizeof(LINK));
    
    for (i = 0; i < g.n; i++) n[i].deg = 0;
    
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
        n[me].deg++;
        n[you].deg++;
    }
    
    for (i = 0; i < g.n; i++) {
        n[i].nb = (NODE **) malloc(n[i].deg * sizeof(NODE *));
        n[i].l = (LINK **) malloc(n[i].deg * sizeof(LINK *));
        n[i].deg = 0;
    }
    rewind(fp);
    
    g.nl = 0;
    while (2 == fscanf(fp, "%u %u\n", &me, &you)) add_link(me, you);
    fclose(fp);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double sir () {
    unsigned int i, now;
    NODE **old = s.i, **new = s.i, *me, *you;
    LINK **lme, **end, **next, *lyou;
    
    // Initialize
    for (i = 0; i < g.n; i++) n[i].t = SUSCEPTIBLE;
    
    // Get a random infection seed, me
    RND_CHK(1);
    me = n + *(g.r++) % g.n;
    
    // Add me as infectious
    *(new++) = me;
    me->t = 0;
    end = s.lchk1;
    for (i = 0; i < me->deg; i++) {
        lyou = me->l[i];
        if (lyou->late == me) {
            you = lyou->late;
            lyou->late = lyou->early;
            lyou->early = you;
        }
        *(end++) = lyou;
    }
    
    // Run the simulation
    for (now = 1; ; now++) {
        // I -> R  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for ( ; old < new; old++) {
            if ((*old)->t + g.delta >= now) break;
            else (*old)->t = RECOVERED;
        }
        
        // if the outbreak has died, break the iterations
        if (new <= old) return (new - s.i) / (double) g.n;
        
        // SI -> II  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // go over the SI, II or IR links
        
        // Chose this time step's link-check list
        if (now & 1) {
            lme = s.lchk1;
            next = s.lchk2;
        } else {
            lme = s.lchk2;
            next = s.lchk1;
        }
        
        RND_CHK(end - lme);
        
        // Go through the link-check list.
        for ( ; lme < end; lme++) {
            you = (*lme)->early;
            
            // Is it an IR link?
            if (you->t != RECOVERED) {
                me = (*lme)->late;
                // Is it an II link?
                if (me->t == SUSCEPTIBLE) {
                    if (*(g.r++) < g.ulambda) {
                        // Infect me
                        *(new++) = me; // Add to the list of infectious nodes
                        me->t = now; // Set the infection times
                        
                        // Checking me's neighbors for more SI links
                        for (i = 0; i < me->deg; i++) {
                            you = me->nb[i]; // (Note that you changes, but that's ok)
                            
                            if (you->t == SUSCEPTIBLE) {
                                lyou = me->l[i];
                                // Add the link to the link-check list
                                *(next++) = lyou;
                                
                                // if needed, change the node order of the link
                                if (lyou->early == you) {
                                    lyou->early = me;
                                    lyou->late = you;
                                }
                            }
                        }
                    } else *(next++) = *lme;
                }
            }
        }
        // Set the end of the next timestep's link list
        end = next;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input

int main (int argc, char *argv[]) {
    unsigned int i;
    double dtmp, s1 = 0.0, s2 = 0.0;
    FILE *fp;
    
    if (argc != 5) {
        fprintf(stderr, "usage: sir_dur3 [rng state file] [nwk file] [lambda] [delta]\n");
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
    if (1 != fread(&g.sfmt, sizeof(sfmt_t), 1, fp)) {
        fprintf(stderr, "error reading rng state\n");
        return 1;
    }
    fclose(fp);
    sfmt_fill_array32(&g.sfmt, g.rnd, NRND);
    g.r = g.rnd;
    // read state or initialize RNG, end
    
    read_network(argv[2]);
    g.lambda = atof(argv[3]);
    // for the transmission probability in units of a 32-bit int
    g.ulambda = (unsigned int) rint(UINT32_MAX * g.lambda);
    g.delta = atoi(argv[4]);
    
    // alloc
    s.i = (NODE **) malloc(g.n * sizeof(NODE *));
    s.lchk1 = (LINK **) malloc(g.nl * sizeof(LINK *));
    s.lchk2 = (LINK **) malloc(g.nl * sizeof(LINK *));
    
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
    if (1 != fwrite(&g.sfmt, sizeof(sfmt_t), 1, fp)) {
        fprintf(stderr, "error writing rng state\n");
        return 1;
    }
    fclose(fp);
    
    // cleaning
    for (i = 0; i < g.n; i++) {free(n[i].nb); free(n[i].l);}
    free(n); free(l); free(g.rnd); free(s.i); free(s.lchk1); free(s.lchk2);
    
    return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
