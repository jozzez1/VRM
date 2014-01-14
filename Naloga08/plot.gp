####
f(x) = k*x + n
fit f(x) "report1.dat" u 1:2 via k, n

plot "report.dat" u 1:2 w l title 'energija', f(x) lt -1 title 'fit'

