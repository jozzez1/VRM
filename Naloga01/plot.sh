# plot.sh
##########

#!/bin/zsh

SAVE=$2

if [ ${SAVE} -eq 0 ]; then

	echo "reset;
	  set ylabel \"x\";
	  set xlabel \"t\";
	  set cblabel \"|u(x)|^2\";
	  set cbrange [0:1];
	  set view map;
	  set palette rgb 34,35,36 maxcolors 10;
	  splot \"$1.txt\" u 1:2:"'($3**2 + $4**2)'" w pm3d title \"|u(x)|^2\" " | gnuplot -p

	  rm $1.txt

else

	echo "reset;
	  set term postscript eps color solid font \"Helvetica,15\";
	  set ylabel \"x\";
	  set xlabel \"t\";
	  set xtics out;
	  set ytics out;
	  set tics nomirror;
	  set cbtics mirror;
	  set cblabel \"|u(x)|^2\";
	  set cbrange [0:1];
	  set view map;
	  set palette rgb 34,35,36 maxcolors 10;
	  set output \"$1.eps\";
	  splot \"$1.txt\" u 1:2:"'($3**2 + $4**2)'" w pm3d title \"|u(x)|^2\" " | gnuplot

	  echo "Converting to a smaller .jpg file"
	  convert -density 220 -quality 100 $1.eps $1.jpg

fi

exit 0

