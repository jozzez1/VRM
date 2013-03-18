# plot3.sh
###########

#!/bin/zsh

DAT="$1"
LAMBDA="$2"
POT="potential-L$LAMBDA.txt"
SAVE="$3"

echo "reset;
set multiplot layout 2,2 rowsfirst;
set autoscale;
set xlabel \"q1\";
set ylabel \"q2\";
set title \"Orbits\";
set view map;
set palette rgbformulae 26,27,28 maxcolor 10;
splot \"$POT\" u 1:2:3 w pm3d title '', \"$DAT.txt\" u 2:3:1 w l lt 5 title '';
set xrange restore;
set yrange restore;
set zrange restore;
set autoscale;
set title 'Energy';
set xlabel 't';
set ylabel 'E';
set grid;
set yrange [0:1];
plot \"$DAT.txt\" u 1:6 w l title 'Total E';
set autoscale;
set title 'index 1';
set xlabel '<p>^2';
set ylabel 't';
plot \"$DAT.txt\" u 1:7 w l title '<p1>^2', \"$DAT.txt\" u 1:8 w l title '<p2>^2';
set autoscale;
set title 'index 2';
set ylabel 'log|E - E0|';
set xlabel 't';
plot \"$DAT.txt\" u 1:"'(log(abs($9)))'" w l title '';
unset multiplot " | gnuplot -p

if [ "$SAVE" == "y" ]; then

	echo -e "Output in $DAT.tex ...\n"
	echo "reset;
	set term epslatex color solid size 14cm, 12cm;
	set output \"$DAT.tex\";
	set multiplot layout 2,2 rowsfirst scale 1.1,1.1;
	set autoscale;
	set xlabel \"q1\";
	set ylabel \"q2\";
	set title \"Orbits\";
	set view map;
	set palette rgbformulae 26,27,28 maxcolor 8;
	splot \"$POT\" u 1:2:3 w pm3d title '', \"$DAT.txt\" u 2:3:1 w l lt 5 title '';
	set xrange restore;
	set yrange restore;
	set zrange restore;
	set autoscale;
	set title 'Energy';
	set xlabel 't';
	set ylabel 'E';
	set grid;
	set yrange [0:1];
	plot \"$DAT.txt\" u 1:6 w l title 'Total E', \"$DAT.txt\" u 1:7 title '<p1>^2', \"$DAT.txt\" u 1:8 title '<p2>^2';
	set autoscale;
	set title 'index 1';
	set xlabel 'p1';
	set ylabel 'q1';
	plot \"$DAT.txt\" u 4:2 w l title '';
	set autoscale;
	set title 'index 2';
	set xlabel 'p2';
	set ylabel 'q2';
	plot \"$DAT.txt\" u 5:3 w l title '';
	unset multiplot " | gnuplot

fi

rm $DAT.txt

exit 0

