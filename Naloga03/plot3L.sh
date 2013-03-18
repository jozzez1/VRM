#!/bin/bash


DAT=$1
SAVE=$2

PIC="pic-$DAT.tex"

echo "reset;
set xlabel 'Lambda';
set ylabel '<p1^2>/<p2^2>'
plot \"$DAT.txt\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot -p

if [ ${SAVE} -eq 1 ]; then

	echo "Saving in $DAT.tex/eps"

	echo "reset;
	set term epslatex color solid size 12cm,8cm;
	set output \"$DAT.tex\";
	set xlabel '$\lambda$';
	set title 'Ekviparticijski test';
	set ylabel '"'$\log\frac{\langle p_1^2\rangle}{/"\langle p_2^2\rangle}$'"
	plot \"$DAT.txt\" u 1:(log("'$2/$3'")) w l title ''" | gnuplot

	epstopdf $DAT.eps
	PIC="joze_zobec_03.tex"

	echo >> $PIC
	echo '\begin{figure}[H]'               >> $PIC
	echo '   \begin{center}'               >> $PIC
	echo '      \input{'"$DAT.tex"'}'      >> $PIC
	echo '   \end{center}'                 >> $PIC
	echo '   \vspace{-20pt}'               >> $PIC
	echo '   \caption{Some caption ...}'   >> $PIC
	echo '   \label{fig:pic-L'"$LAMBDA"'}' >> $PIC
	echo '   \vspace{-10pt}'               >> $PIC
	echo '\end{figure}'                    >> $PIC
	echo >> $PIC

fi

exit 0

