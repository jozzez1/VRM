function U = potentialU (u, rmax, N)

	h = rmax / N;

	r = linspace (1e-5, rmax, N);
	ur = - (u.^2) ./ r;
	L = mLaplace (rmax, N);
	
	[nr, nc] = size (ur);
	if (nc != 1)
		ur = ur.';
	endif
	U = L \ ur;

	% we still have to add the linear part
	k = (1 - U(N))/rmax

	U = U + k*r.';

	return;
endfunction
