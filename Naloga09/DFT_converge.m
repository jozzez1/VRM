function [E2, u2, U2] = DFT_converge (E, u, U, rmin, rmax)
	E1 = E;
	u1 = u;
	U1 = U;

	[E2, u2, U2] = DFT_step (E1, u1, U1, rmin, rmax);
	while (sumsq(u2 - u1) + sumsq(U2 - U1) > 1e-6)
		E1 = E2;
		u1 = u2;
		U1 = U2;

		[E2, u2, U2] = DFT_step (E1, u1, U1, rmin, rmax);
	endwhile

	return;
endfunction
