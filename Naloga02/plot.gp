reset

set term epslatex color solid size 14cm,7cm
set output "spectre.tex"

set xlabel '$N$'
set ylabel '$E_N$

set key top left

set log xy

f(x) = x-0.5
plot "Energies-N2000-L0.txt" u 1:2 w l lt 1 title '$\lambda = 0.01$', \
f(x) lt -1 w l title '$\lambda = 0$ pravi', \
"Energies-N2000-L1.txt" u 1:2 w l lt 2 title '$\lambda = 0$', \
"Energies-N2000-L100.txt" u 1:2 w l lt 3 title '$\lambda = 1$'
