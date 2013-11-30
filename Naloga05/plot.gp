# temperature profiles
######################

#1.1
reset
set term epslatex color solid
set output "temp-N20.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle T_i \rangle$'
set key top left
set title 'Temperaturni profil za $N = 20$'

plot "T-N20-L0.txt" w lp title '$\lambda = 0$', \
	"T-N20-L200.txt" w lp title '$\lambda = 2$', \
	"T-N20-L400.txt" w lp title '$\lambda = 4$'

unset output

#2.1
reset
set term epslatex color solid
set output "temp-N40.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle T_i\rangle $'
set key top left
set title 'Temperaturni profil za $N = 40$'

plot "T-N40-L0.txt" w lp title '$\lambda = 0$', \
	"T-N40-L200.txt" w lp title '$\lambda = 2$', \
	"T-N40-L400.txt" w lp title '$\lambda = 4$'

unset output

#3.1
reset
set term epslatex color solid
set output "temp-N60.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle T_i\rangle $'
set key top left
set title 'Temperaturni profil za $N = 60$'

plot "T-N60-L0.txt" w lp title '$\lambda = 0$', \
	"T-N60-L200.txt" w lp title '$\lambda = 2$', \
	"T-N60-L400.txt" w lp title '$\lambda = 4$'

unset output

#4.1
reset
set term epslatex color solid
set output "temp-N80.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle T_i\rangle $'
set key top left
set title 'Temperaturni profil za $N = 80$'

plot "T-N80-L0.txt" w lp title '$\lambda = 0$', \
	"T-N80-L200.txt" w lp title '$\lambda = 2$', \
	"T-N80-L400.txt" w lp title '$\lambda = 4$'

unset output

# current profiles
######################

#1.2
reset
set term epslatex color solid
set output "curr-N20.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle J_i\rangle $'
set title 'Energijski tok za $N = 20$'
set key out

plot "J-N20-L0.txt" w lp title '$\lambda = 0$', \
	"J-N20-L200.txt" w lp title '$\lambda = 2$', \
	"J-N20-L400.txt" w lp title '$\lambda = 4$'

#2.2
reset
set term epslatex color solid
set output "curr-N40.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle J_i\rangle $'
set title 'Energijski tok za $N = 40$'
set key out

plot "J-N40-L0.txt" w lp title '$\lambda = 0$', \
	"J-N40-L200.txt" w lp title '$\lambda = 2$', \
	"J-N40-L400.txt" w lp title '$\lambda = 4$'

#3.2
reset
set term epslatex color solid
set output "curr-N60.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle J_i\rangle $'
set title 'Energijski tok za $N = 60$'
set key out

plot "J-N60-L0.txt" w lp title '$\lambda = 0$', \
	"J-N60-L200.txt" w lp title '$\lambda = 2$', \
	"J-N60-L400.txt" w lp title '$\lambda = 4$'

#4.2
reset
set term epslatex color solid
set output "curr-N80.tex"
set xlabel 'indeks \v clenka verige -- $i$'
set ylabel '$\langle J_i\rangle $'
set title 'Energijski tok za $N = 40$'
set key out

plot "J-N80-L0.txt" w lp title '$\lambda = 0$', \
	"J-N80-L200.txt" w lp title '$\lambda = 2$', \
	"J-N80-L400.txt" w lp title '$\lambda = 4$'

# current average fit
reset
set term epslatex color solid
set output "avg-curr.tex"
set xlabel 'dol\v zina verige -- $N$'
set ylabel '$\langle J \rangle$'
set title 'Povpre\v cni $\langle J_j\rangle $ za $\lambda = 4$'

f (x) = k / x;
fit f(x) "avg-J-L4.txt" u 1:2:3 via k
plot "avg-J-L4.txt" u 1:2:3 w yerror title 'podatki', \
	f(x) title 'fit' lt -1
unset output

set output "avg-curr2.tex"
set key top left
plot "avg-J-L4.txt" u 1:2 w lp title 'meritve',\
	f(x) title 'fit' lt -1
unset output
