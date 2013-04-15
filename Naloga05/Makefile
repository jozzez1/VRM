CC=clang
CFLAGS=-O3 -march=corei7 -lm -Wall -pedantic-errors -lc
PRGNAM=nal5
SOURCE=prg5.c

all:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS)

dbg:
	$(CC) -o $(PRGNAM) $(SOURCE) $(CFLAGS) -g

clean:
	rm -rf T* J* $(PRGNAM)*
