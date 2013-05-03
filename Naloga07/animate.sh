#!/bin/sh

dir=$1
base=$2
#range=`echo "scale = 1; $3-1" | bc`
range=$3
temperature=$4

font=/usr/local/lib/X11/fonts/ChromeOS/Arimo-Regular.ttf

echo "set terminal jpeg enhanced font \"$font,12\" size 550,400;
  unset surf;
  set pm3d;
  set view map;
  set palette rgbformulae 21,22,23 maxcolor 10;
  set title \"Temperature = $temperature\";
  set output \"$dir/$base.jpg\";
  set xrange [0:$range];
  set yrange [0:$range];
  set cbrange [-0.8:0.8];
  plot \"$dir/$base.txt\" matrix w image title ''" | gnuplot

exit 0

