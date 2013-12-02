CC=clang
PRGNAM=nal7
SOURCE=prg7.c
CPUTYPE=corei7
COPTFLAG?=-O3

CFLAGS=$(COPTFLAG) -march=$(CPUTYPE) -pipe\
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
	rm -rf video-* N*-M*.txt
