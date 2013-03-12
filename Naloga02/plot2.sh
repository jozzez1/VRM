#!/bin/zsh

MODE=$1
DATA=$2
SAVE=$3

if [ ${MODE} -eq 0 ]; then

	if [ "$SAVE" != "1" ]; then
		echo "reset;
		set ylabel 'h';
		set xlabel 'N';
		set y2label 'good';
		set tics nomirror;
		set autoscale y2;
		set y2tics;
		f(x) = a/(x**b);
		fit f(x) \"h-of-N.txt\" via a, b;
		plot \"h-of-N.txt\" u 1:2 title '', f(x) lt -1 title 'fit', \
		\"h-of-N.txt\" u 1:3 lt 3 title 'good' axes x1y2" | gnuplot -p
	
	else
		echo "reset;
		set ylabel 'h';
		set xlabel 'N';
		set y2label 'good';
		set tics nomirror;
		set autoscale y2;
		set y2tics;
		f(x) = a/(x**b);
		fit f(x) \"h-of-N.txt\" via a, b;
		plot \"h-of-N.txt\" u 1:2 title '', f(x) lt -1 title 'fit', \
		\"h-of-N.txt\" u 1:3 lt 3 title 'good' axes x1y2" | gnuplot -p

		echo "reset;
		set term epslatex color solid size 8cm,4cm;
		set output \"calibrate.tex\";
		set ylabel 'h';
		set xlabel 'N';
		set y2label 'good';
		set tics nomirror;
		set autoscale y2;
		set y2tics;
		f(x) = a/(x**b);
		fit f(x) \"h-of-N.txt\" via a, b;
		plot \"h-of-N.txt\" u 1:2 title '', f(x) lt -1 title 'fit', \
		\"h-of-N.txt\" u 1:3 lt 3 title 'good' axes x1y2" | gnuplot -p

	fi

else

	echo "reset;
	set key top left;
	set ylabel 'r';
	set xlabel 'N';
	f(x) = a*(x**b);
	fit f(x) \"ratio-L1-fit.txt\" u 1:2 via a, b;
	plot \"$DATA\" u 1:2 title '', f(x) lt -1 title 'fit'" | gnuplot -p

	if [ "$SAVE" != "1" ]; then
		echo "reset;
		set key top left;
		set term epslate color solid size 8cm,4cm;
		set output \"Lambda.tex\";
		set xlabel 'N';
		set ylabel 'r';
		f(x) = a*(x**b);
		fit f(x) \"ratio-L1-fit.txt\" u 1:2 via a, b;
		plot \"$DATA\" u 1:2 title '', f(x) lt -1 title 'fit'" | gnuplot -p
	fi

fi

exit 0

