# FreeBSD 9.1 compliant Makefile
##################################
CC=gcc48
CFLAGS=-Wall -lm -O3 -march=native -lgsl -latlas -lalapack -lptcblas -lf77blas -lc
PRGNAM=nal2
SOURCE=prg2.c

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g
