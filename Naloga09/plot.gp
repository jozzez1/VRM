reset

set term epslatex color solid

set output "energio.tex"
set title '$u_N(E)$'
set xlabel 'energija -- $E$'
set ylabel '$u(r_\text{max})$'
set grid

plot "energio.txt" u 1:2 w l title '$u$'
unset output

set output "vodik.tex"
set title 'Rezultati za vodikov atom'
set xlabel 'Radialna koordinata -- $r$'
set ylabel 'Amplituda'

plot "vodik.txt" u 1:2 w l title '$u(r)$', \
"vodik.txt" u 1:3 w l lt 3 title '$U(r)$'
unset output

set output "helij.tex"
set title 'Rezulztati za helijev atom'

plot "helij.txt" u 1:2 w l title '$u(r)$', \
"helij.txt" u 1:3 w l lt 3 title '$U(r)$'
unset output


set output "uerr.tex"
set ylabel '$\log_{10}|U - U_\text{ex}|$'
set title 'Ujemanje za analiti\v cno vrednostjo'

plot "vodik.txt" u 1:(log10(abs($3 - $4))) w l title 'Logaritem napake'
unset output

