reset

set term epslatex color solid
set output "FreeEn.tex"
set xlabel '$\beta$'
set ylabel '$F(\beta)$'

plot "F-G5-N20.txt" u 1:3 w l title '$N = 20$', \
"F-G5-N16.txt" u 1:3 w l title '$N = 16$', \
"F-G5-N12.txt" u 1:3 w l title '$N = 12$', \
"F-G5-N8.txt" u 1:3 w l title '$N = 8$', \
"F-G5-N6.txt" u 1:3 w l title '$N = 6$', \
"F-G5-N2.txt" u 1:3 w l title '$N = 2$'

set output "Energy.tex"
set ylabel '$\langle H \rangle_\beta$'

plot "E-G5-N20.txt" u 1:4 w l title '$N = 20$', \
"E-G5-N16.txt" u 1:4 w l title '$N = 16$', \
"E-G5-N12.txt" u 1:4 w l title '$N = 12$', \
"E-G5-N8.txt" u 1:4 w l title '$N = 8$', \
"E-G5-N6.txt" u 1:4 w l title '$N = 6$', \
"E-G5-N2.txt" u 1:4 w l title '$N = 2$'


