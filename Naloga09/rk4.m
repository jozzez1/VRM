function unew = rk4 (V, E, rmax, N)
	% initial conditions
	u{1} = [-1/V(1), 0.1];

	% we set some parameters
	f = 2 .* (V .- E);
	h = (rmax - 1e-5)/N;

	% here we use the actual RK4 method
	for k = 1:N-1
		% full-step A
		A1 = [0, 1;f(k), 0];
		
		% half-step A
		A2 = [0, 1; 0.5*(f(k) + f(k+1)), 0];

		% two-step A
		A3 = [0, 1; f(k+1), 0];

		k1 = A1 * u{k}.';
		k2 = A2 * (u{k}.' + 0.5*h .* k1);
		k3 = A2 * (u{k}.' + 0.5*h .* k2);
		k4 = A3 * (u{k}.' + h .* k3);

		u{k+1} = u{k} + (h/6) .* (k1.' + 2*k2.' + 2*k3.' + k4.');
	end

	% all we need to do now is to take the
	% u(r) and discard the gradient, u'(r)
	for k = 1:N
		unew (k) = u{k}(1);
	end

	% u(r) must be normalized the integral sense
	unew = sqrt(N) * unew ./ norm(unew);

	return;
endfunction
