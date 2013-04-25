#!/bin/sh

dir=$1
base=$2
range=`echo "scale = 1; $3-1" | bc`
temperature=$4

font=/usr/local/lib/X11/fonts/ChromeOS/Arimo-Regular.ttf

echo "set terminal jpeg enhanced font \"$font,12\" size 550,400;
  unset surf;
  set view map;
  set palette rgbformulae 10,10,10 maxcolor 2;
  set title \"Temperature = $temperature\";
  set output \"$dir/$base.jpg\";
  set xrange [-1:$range];
  set yrange [-1:$range];
  set cbrange [-1.1:1.1];
  plot \"$dir/$base.txt\" matrix with image title ''" | gnuplot

exit 0

