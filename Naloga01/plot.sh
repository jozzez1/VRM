# plot.sh
##########

#!/bin/zsh

SAVE=$2

if [ ${SAVE} -eq 0 ]; then

	echo "reset;
	  set ylabel \"x\";
	  set xlabel \"t\";
	  set cblabel \"|u(x)|^2\";
	  set view map;
	  splot \"$1.txt\" u 1:2:"'($3**2 + $4**2)'" w pm3d title \"|u(x)|^2\" " | gnuplot -p

	  rm $1.txt

else

	echo "reset;
	  set term epslatex color solid;
	  set ylabel \"x\";
	  set xlabel \"t\";
	  set cblabel \"|u(x)|^2\";
	  set view map;
	  set output \"$1.eps\";
	  splot \"$1.txt\" u 1:2:"'($3**2 + $4**2)'" w pm3d title \"|u(x)|^2\" " | gnuplot

fi

exit 0

