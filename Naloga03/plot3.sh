# plot3.sh
###########

#!/bin/zsh

DAT="$1"
LAMBDA="$2"
POT="potential-L$LAMBDA.txt"
SAVE="$3"

PIC="pic-$DAT.tex"

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
	set multiplot layout 2,2 rowsfirst;
	set autoscale;
	set xlabel '$"'q_1'"$';
	set ylabel '$"'q_2'"$';
	set title \"Orbita\";
	set view map;
	set palette rgbformulae 26,27,28 maxcolor 10;
	splot \"$POT\" u 1:2:3 w pm3d title '', \"$DAT.txt\" u 2:3:1 w l lt 5 title '';
	set xrange restore;
	set yrange restore;
	set zrange restore;
	set autoscale;
	set title 'Energija';
	set xlabel '$"'t'"$';
	set ylabel '"'$E'"';
	set grid;
	set xtics 200;
	set yrange [0:1];
	plot \"$DAT.txt\" u 1:6 w l title 'celotna E';
	set autoscale;
	set title 'Povprečja';
	set xlabel '$\langle p^2 "'\r'"angle$';
	set ylabel '"'$t'"$';
	set xtics 200;
	plot \"$DAT.txt\" u 1:7 w l title '$\langle p_1^2 "'\r'"angle$', \"$DAT.txt\" u 1:8 w l title '$\langle p_2^2 "'\r'"angle$';
	set autoscale;
	set xtics 200;
	set title 'Natančnost energije';
	set ylabel '$\log|E - E_0|$';
	set xlabel '"'$t'"$';
	plot \"$DAT.txt\" u 1:"'(log(abs($9)))'" w l title '';
	unset multiplot " | gnuplot

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

	rm $DAT.txt

fi

exit 0

