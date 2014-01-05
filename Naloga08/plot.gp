# gnuplot file ...
#1st let's get the compact

set xlabel 'dol\v zina verige -- $n$'
set ylabel 'Entropija prepletenosti -- $S$'

set key top left
set term epslatex color solid

set output "1-kompakten.tex"
set title 'Kompaktna izbira biparticije'
plot "1-non-and-compact.txt" u 1:4 w lp title 'odprt', \
"1-non-and-compact.txt" u 1:5 w lp lt 3 title 'periodi\v cen'
unset output

set output  "1-nekompakten.tex"
set title 'Nekompaktna izbira biparticije'
plot "1-non-and-compact.txt" u 1:2 w lp title 'odprt', \
"1-non-and-compact.txt" u 1:3 w lp lt 3 title 'periodi\v cen'
unset output

# now let's plot the variations over 'a'

set xrange [0:1]
set xlabel 'Dele\v z spinov v levi particiji -- $|A|/n$'
set title 'Entropija v odvisnosti od velikosti leve particije'
set key out top right

set output "1-periodicni-aji.tex"
plot "1-periodic.txt" u 1:2 w lp title '$n = 2$', \
"1-periodic.txt" u 3:4 w lp title '$n = 4$', \
"1-periodic.txt" u 5:6 w lp title '$n = 6$', \
"1-periodic.txt" u 7:8 w lp title '$n = 8$', \
"1-periodic.txt" u 9:10 w lp title '$n = 10$'
unset output

set output "1-neperiodicni-aji.tex"
plot "1-nonperiod.txt" u 1:2 w lp title '$n = 2$', \
"1-nonperiod.txt" u 3:4 w lp title '$n = 4$', \
"1-nonperiod.txt" u 5:6 w lp title '$n = 6$', \
"1-nonperiod.txt" u 7:8 w lp title '$n = 8$', \
"1-nonperiod.txt" u 9:10 w lp title '$n = 10$'
unset output

set output "1-ni-parabola-a.tex"
set key inside top right
f(x) = a*x*(x-1)
g(x) = b*(x-0.5)**4 + c
fit f(x) "1-periodic.txt" u 9:10 via a
fit g(x) "1-periodic.txt" u 9:10 via b,c
plot "1-periodic.txt" u 9:10 w lp title 'n = 12', f(x) lt -1 title 'fit $\propto x^2$', g(x) lt 3 title 'fit $\propto x^4$'
unset output

