function E0 = energy (E, u, U, rmin, rmax, N)
	N = length (u);
	r = linspace (rmin, rmax, N);
	h = (rmax - rmin)/N;

	Vhf = 2 * (U ./ r);
	Vxc = (-1)*cbrt(3/(2*pi^2)) * cbrt((u.^2)./(r.^2));


	E0 = 2*E - h*(sum(Vhf .* (u.^2)) + 0.5 * sum(Vxc .* (u.^2)));

	return;
endfunction
