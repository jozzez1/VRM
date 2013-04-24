CC=gcc48
PRGNAM=nal6
SOURCE=prg6.c
CPUTYPE=native

CFLAGS=-O3 -march=$(CPUTYPE) -pipe\
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
	rm -rf J* ISING*
