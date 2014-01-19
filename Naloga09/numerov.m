function u = numerov (V, E, rmin, rmax, N)
	% initial conditions
	u(1) = 0;
	u(2) = 1;

	% ok ... we set some parameters
	f = 2 * (E .- V);
	h = (rmax - rmin)/N;

	% and here we do the iteration ...
	for k = 2:N-1
		u(k+1) = (2 - (5/6)*h^2*f(k))*u(k);
		u(k+1) -= (1 + (f(k-1)*h^2)/12)*u(k-1);
		u(k+1) /= 1 + (f(k+1)*h^2)/12;
	end

	u = u ./ norm(u);
	u *= sqrt(N/(rmax - rmin));

	return;
endfunction
