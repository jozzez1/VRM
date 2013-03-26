#!/bin/zsh

DAT=$1
SAVE=$2
MODE=$3
NUM=$4

# write stuff in my report ...
function write ()
{
	dat=$1
	report=$2

	# we don't want to change this all the time ...
	sed -i '/\end{document}/d' $report

	echo -E '\begin{figure}[H]'               >> $report
	echo -E '   \begin{center}'               >> $report
	echo -E '      \input{'"$dat.tex"'}'      >> $report
	echo -E '   \end{center}'                 >> $report
	echo -E '   \vspace{-20pt}'               >> $report
	echo -E '   \caption{Some caption ...}'   >> $report
	echo -E '   \label{fig:pic}'              >> $report
	echo -E '   \vspace{-10pt}'               >> $report
	echo -E '\end{figure}'                    >> $report
	echo                                      >> $report
	echo -E '\end{document}'                  >> $report
	echo >> $report
}

function plot ()
{
	save=$1
	mode=$2
	dat=$3
	num=$4

	temp=scriptP.gp

	rm -rf $temp
	touch $temp

	set pflag
	set xaxis
	set yaxis
	set y2axis

	if [ ${mode} -eq 0 ]
	then
		xaxis='$\beta$'
		yaxis='$F$'
		y2axis='$Z$'
	elif [ ${mode} -eq 1 ]
	then
		xaxis='$\beta$'
		yaxis='$E$'
	elif [ ${mode} -eq 2 ]
	then
		xaxis='$t$'
		yaxis='$C$'
	elif [ ${mode} -eq 3 ]
	then
		xaxis='$t$'
		yaxis='$J$'
	fi

	# we prepare the terrain
	echo >> $temp
	echo 'reset' >> $temp

	if [ ${save} -eq 1 ]
	then
		echo 'set term epslatex color solid' >> $temp
		echo "set output \"$dat.tex\""       >> $temp

	else
		pflag="-p"
	fi

	if [ ${mode} -eq 0 ]; then

		echo 'set tics nomirror'    >> $temp
		echo 'set autoscale y2'     >> $temp
		echo 'set y2tics'           >> $temp
		echo 'set log y2'           >> $temp
		echo "set y2labe '$y2axis'" >> $temp
	fi

	echo -E "set xlabel '$xaxis'"  >> $temp
	echo -E "set ylabel '$yaxis'"  >> $temp
	echo -E "set label '"'$N = '"$num"'$'"' at screen 0.72,0.3" >> $temp

	# here comes the plotting part
	echo >> $temp
	if [ ${mode} -eq 0 ]; then
		echo -E "plot \"$dat.txt\" u 1:3 w l title '$yaxis' axes x1y1, \\" >> $temp
		echo -E "\"$dat.txt\" u 1:2 w l lt 3 title '$y2axis' axes x1y2"    >> $temp

	else
		echo -E "plot \"$dat.txt\" u 1:4 w l title '$yaxis'"               >> $temp
	fi

	# we finish with the file ...
	echo >> $temp

	# ... so we start plotting ... ;)
	cat $temp | gnuplot $pflag
}

if [ ${MODE} -ne 0 ]
then
	plot 0 0 $DAT $NUM
fi

plot 0 $MODE $DAT $NUM

if [ ${SAVE} -eq 1 ]; then
	plot $SAVE $MODE $DAT $NUM
	epstopdf $DAT.eps
	write $DAT joze_zobec_04.tex

	if [ ${MODE} -eq 0 ]; then
		plot $SAVE 0 $MODE $DAT $NUM
		epstopdf $DAT.eps
		write $DAT joze_zobec_04.tex
	fi
fi

exit 0

