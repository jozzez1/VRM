reset
set term epslatex color solid size 10cm,11cm
set tics out
set tics nomirror
set view map
set yrange [] reverse
set palette rgbformulae 21,22,23 maxcolors 10

# dolzina 10
set output "korelacija-10-0.tex"
splot "korelacija-10-0.txt" matrix w image title ''
unset output

# dolzina 12
set output "korelacija-12-0.tex"
splot "korelacija-12-0.txt" matrix w image title ''
unset output
