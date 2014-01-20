rmax = 15;
rmin = 1e-20;
N = 8000;

E = linspace (-15, 5, 1000);
r = linspace (rmin, rmax, N);

V = -1 ./ r;

for k = 1:1000
	y(k) = numerov (V, E(k), rmin, rmax, N)(N);
end

R (:,1) = E.';
R (:,2) = y.';

save -ascii energio.txt R
