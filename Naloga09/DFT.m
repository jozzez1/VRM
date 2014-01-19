function [E, u, U] = DFT (rmin, rmax, N)
	% we initialize our starting vectors
	[E, u, U] = start (-0.5, rmin, rmax, N);

	% now we iterate it
	[E, u, U] = DFT_converge (E, u, U, rmin, rmax);

	% well, now we could integrate it ... should we?

	return;
	
endfunction
