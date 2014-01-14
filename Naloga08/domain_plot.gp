reset
set view map
set tics nomirror
set tics out
set title 'Razvoj domenske stene, $n = 12$, $M = 0$'
set xlabel 'indeks spina v verigi -- $j$'
set cblabel '$\langle\sigma^z_j\rangle$'
set ylabel 'time -- $t$'
#set cbrange [-1:1]
set term epslatex color solid size 10cm,14cm
set palette rgbformulae 30,31,32 maxcolors 40

set ytics ("0" 0, "2" 100, "4" 200, "6" 300, "8" 400, "10" 500)
set xtics 1
set xrange [-0.5:11.5]

set cbtics mirror
set cbtics in

set output "profil-12-32.tex"
set title 'Razvoj domenske stene, $n = 12$, $M = 32$'
splot "domain-profile-12-32.txt" matrix w image title ''
unset output

set cbtics 0.25
set output "profil-12-0.tex"
splot "domain_profile.txt" matrix w image title ''
unset output

reset
set term epslatex color solid size 14cm,10cm
set title 'Posamezni spini'
set xlabel '\v cas -- $t$'
set ylabel '$\langle\sigma_j^z\rangle$'
set tics out
set tics nomirror
set grid
set key b r
set output "spini.tex"
plot "domain.txt" u 1:2 w l title '$j = 0$', \
"domain.txt" u 1:3 w l title '$j = 1$', \
"domain.txt" u 1:4 w l title '$j = 2$', \
"domain.txt" u 1:5 w l title '$j = 3$', \
"domain.txt" u 1:6 w l title '$j = 4$', \
"domain.txt" u 1:7 w l title '$j = 5$' lt -1
unset output

