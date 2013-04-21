#!/usr/local/bin/zsh

echo -e "reset;
   unset surf;
   set view map;
   set palette rgbformulae 10,10,10 maxcolor 2;
   splot \"$1.txt\" matrix with image;" | gnuplot -p
