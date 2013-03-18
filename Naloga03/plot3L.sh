#!/bin/zsh

DAT=$1
TAU=$2
TTT=$3
SAVE=$4

echo "reset;
set xlabel 'Lambda';
set ylabel '<p1^2>/<p2^2>'
plot \"$DAT\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot -p

if [ "$SAVE" == "y" ]; then
	echo "reset;
	set term epslatex color solid size 12cm,8cm;
	set output 'scan-$TAU-$TTT.tex';
	set xlabel '$\lambda$';
	set ylabel '$"'\r'"angle p_1^2\langle/"'\r'"anglep_2^2\langle'
	plot \"$DAT\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot -p
fi

exit 0

