#!/bin/zsh


DAT=$1
SAVE=$2

echo "reset;
set xlabel 'Lambda';
set ylabel '<p1^2>/<p2^2>'
plot \"$DAT.txt\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot -p

if [ ${SAVE} -eq 1 ]; then

	echo -e "Saving in $DAT.tex/eps\n"
	echo "reset;
	set term epslatex color solid size 12cm,8cm;
	set output \"$DAT.tex\";
	set xlabel '$\lambda$';
	set ylabel '$"'\r'"angle p_1^2\langle/"'\r'"anglep_2^2\langle'
	plot \"$DAT.txt\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot -p

fi

exit 0

