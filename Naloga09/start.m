function [E, u, U] = start (E0, rmin, rmax, N)
	r = linspace (rmin, rmax, N);
	V = -1./r;

	[E, u] = ground (V, E0, rmin, rmax, N);
	U = potentialU (u, rmin, rmax, N);

	return;
endfunction
