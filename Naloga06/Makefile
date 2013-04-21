CC=clang
PRGNAM=nal6
SOURCE=prg6.c
CPUTYPE=corei7-avx

CFLAGS=-O3 -march=$(CPUTYPE) -pipe\
       -Wall -pedantic-errors\
       -lm -lc -lgsl -lgslcblas\
       -I/usr/local/include\
       -L/usr/local/lib\
       -std=c99

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g

clean:
	rm -rf c0* c1*
