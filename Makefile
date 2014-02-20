SRC = .
SFMT_SRC = .
CFLAGS = -W -Wall -DSFMT_MEXP=19937 -DSFMT_HAVE_SSE2=1 -Ofast -msse2 -march=native
LDFLAGS = 
CC = gcc

OBJ1 = o/sir_dur.o o/SFMT.o
OBJ2 = o/sir_dur2.o o/SFMT.o
OBJ3 = o/sir_dur3.o o/SFMT.o

all : sir_dur sir_dur2 sir_dur3

sir_dur: $(OBJ1)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

sir_dur2: $(OBJ2)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

sir_dur3: $(OBJ3)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

o/SFMT.o : $(SFMT_SRC)/SFMT.c $(SFMT_SRC)/SFMT.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SFMT_SRC)/SFMT.c -o $@

o/sir_dur.o : $(SRC)/sir_dur.c $(SRC)/sir_dur.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/sir_dur.c -o $@

o/sir_dur2.o : $(SRC)/sir_dur2.c $(SRC)/sir_dur2.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/sir_dur2.c -o $@

o/sir_dur3.o : $(SRC)/sir_dur3.c $(SRC)/sir_dur3.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/sir_dur3.c -o $@
