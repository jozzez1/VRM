function [Enew, unew, Unew] = DFT_step (E, u, U, rmin, rmax)

	N = length (u);
	r = linspace (rmin, rmax, N);

	Vhf = 2*U./r;
	Vxc = (-1)*cbrt(3/(2*pi^2)) * cbrt((u.^2)./(r.^2));

	V = (-2)./r + Vhf + Vxc;

	[Enew, unew] = ground (V, E, rmin, rmax, N);
	Unew = potentialU (unew, rmin, rmax, N);

	return;

endfunction
