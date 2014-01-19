function [E0, E, u, U] = DFT (rmin, rmax, N, prec)
	% we initialize our starting vectors
	[E, u, U] = start (-0.5, rmin, rmax, N);

	% now we iterate it
	[E, u, U] = DFT_converge (E, u, U, rmin, rmax, prec);

	% well, now we could integrate it ... should we?
	E0 = energy (E, u, U, rmin, rmax, N);

	return;
	
endfunction
