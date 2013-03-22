#!/bin/bash

DAT=$1
SAVE=$2
MODE=$3
NUM=$4

if [ ${MODE} -eq 0 ]; then

	echo "reset;
	set tics nomirror;
	set autoscale y2;
	set key top left;
	set label \"N = $NUM\" at screen 0.72,0.3;
	set y2tics;
	set log y2;
	set ylabel 'F';
	set y2label 'Z';
	plot \"$DAT.txt\" u 1:3 w l title 'F' axes x1y1;
	replot \"$DAT.txt\" u 1:2 w l lt 3 title 'Z' axes x1y2" | gnuplot -p
	
	if [ ${SAVE} -eq 1 ]; then
	
		echo "reset;
		set term epslatex color solid;
		set output \"$DAT.tex\";
		set tics nomirror;
		set autoscale y2;
		set key top left;
		set y2tics;
		set log y2;
		set tics out;
		set label \"N = $NUM\" at screen 0.72,0.3;
		set title 'Prosta energija';
		set xlabel '"'$\beta$'"';
		set ylabel '"'$F$'"';
		set y2label '"'$Z$'"';
		plot \"$DAT.txt\" u 1:3 w l title '"'$F(\beta)$'"' axes x1y1;
		replot \"$DAT.txt\" 1:2 w l lt 3 title '"'$Z(\beta)$'"' axes x1y2 "| gnuplot
	
		epstopdf $DAT.eps
		PIC="joze_zobec_04.tex"
	
		echo >> $PIC
		echo '\begin{figure}[H]'               >> $PIC
		echo '   \begin{center}'               >> $PIC
		echo '      \input{'"$DAT.tex"'}'      >> $PIC
		echo '   \end{center}'                 >> $PIC
		echo '   \vspace{-20pt}'               >> $PIC
		echo '   \caption{Some caption ...}'   >> $PIC
		echo '   \label{fig:pic}'              >> $PIC
		echo '   \vspace{-10pt}'               >> $PIC
		echo '\end{figure}'                    >> $PIC
		echo >> $PIC

	fi

elif [ ${MODE} -eq 1 ]; then

	echo "reset;
	set tics nomirror;
	set autoscale y2;
	set key top left;
	set y2tics;
	set ylabel 'Re (E)'
	set y2label 'Im (E)'
	set label \"N = $NUM\" at screen 0.72,0.3;
	set log y2;
	plot \"$DAT.txt\" u 1:4 w l title 'Re (E)' axes x1y1;
	replot \"$DAT.txt\" u 1:5 w l lt 3 title 'Im (E)' axes x1y2" | gnuplot -p

elif [ ${MODE} -gt 1 ]; then

	echo "reset;
	set ylabel 'C';
	set xlabel 't';
	set label \"N = $NUM\" at screen 0.72,0.3;
	plot \"$DAT.txt\" u 1:2 w l title 'Re (C)';
	replot \"$DAT.txt\" u 1:3 w l lt 3 title 'Im (C)'" | gnuplot -p

fi

#rm -rf $DAT.txt

exit 0

