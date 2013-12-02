CC=cc
PRGNAM=nal5
SOURCE=prg5.c
CPUTYPE=native

CFLAGS=-O3 -march=$(CPUTYPE) -pipe\
       -pedantic-errors\
       -lm -lc -lgsl -lgslcblas\
       -I/usr/local/include\
       -L/usr/local/lib\
       -std=c99

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g

clean:
	rm -rf T* J* H*
