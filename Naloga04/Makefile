CC=cc
CFLAGS=-Wall -lm -lgsl -lgslcblas -march=native -O3
PRGNAM=nal4
SOURCE=prg4.c

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g
