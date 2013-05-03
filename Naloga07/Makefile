CC=gcc
PRGNAM=nal7
SOURCE=prg7.c
CPUTYPE=native
COPTFLAG?=-O3

CFLAGS=$(COPTFLAG) -march=$(CPUTYPE) -pipe\
       -Wall -pedantic-errors\
       -lm -lc -lgsl -lgslcblas -fopenmp\
       -I/usr/local/include\
       -L/usr/local/lib\
       -std=c99

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g

clean:
	rm -rf animate* HARMONIC*
