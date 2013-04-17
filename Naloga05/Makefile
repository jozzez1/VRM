CC=clang
PRGNAM=nal5a
SOURCE=prg5.c

CFLAGS=-O3 -march=corei7-avx -Wall -pedantic-errors -lm -lc\
       -lgsl -lgslcblas\
       -I/usr/local/include\
       -L/usr/local/lib

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g

clean:
	rm -rf T* J* $(PRGNAM).core
