#!/bin/zsh

MODE=$1
DATA=$2
SAVE=$3

if [ ${MODE} -eq 0 ]; then

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
	set ylabel 'r';
	set xlabel 'N';
	f(x) = a/(x**b) - c*x;
	plot \"$DATA\" u 1:2 title ''" | gnuplot -p

	if [ "$SAVE" != "1" ]; then
		rm $DATA
	fi

fi
