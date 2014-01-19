function [E0, u] = ground (V, rmax, N)
	E1 = -10;
	E2 = 0.02;

	u1 = rk4 (V, E1, rmax, N);
	u2 = rk4 (V, E2, rmax, N);
	while (1e-8 < abs(E1 - E2))

		E3 = (E1 + E2)/2;
		u3 = rk4 (V, E3, rmax, N);

		if u3(N) * u1(N) < 0
			E2 = E3;
			u2 = u3;
		elseif u3(N)*u2(N) < 0
			E1 = E3;
			u1 = u3;
		else
			printf ("Error!\n");
			break;
		endif
	end

	E0 = E1;
	u  = u1;

	return;
endfunction
