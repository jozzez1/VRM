E0 = -0.5;
rmin = 1e-20
rmax = 15
N = 12000

[E, u, U] = start (E0, rmin, rmax, N);

printf ("E = %f\n", E);

r = linspace (rmin, rmax, N);
f = -(r .+ 1) .* exp(-2*r) .+ 1;

R (:,1) = r.';
R (:,2) = u.';
R (:,3) = U.';
R (:,4) = f.';

save -ascii vodik.txt R

