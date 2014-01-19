function [E0, u] = ground (V, E, rmin, rmax, N)
	E1 = E-2;
	E2 = E+2;

	u1 = numerov (V, E1, rmin, rmax, N);
	u2 = numerov (V, E2, rmin, rmax, N);

	counter = 0;
	while (1e-10 < abs(E1 - E2))

		E3 = (E1 + E2)/2;
		u3 = numerov (V, E3, rmin, rmax, N);

		if u3(N) * u1(N) < 0
			E2 = E3;
			u2 = u3;
		elseif u3(N)*u2(N) < 0
			E1 = E3;
			u1 = u3;
		elseif u1(N)*u2(N) > 0
			while u1(N)*u2(N) > 0
				E1 += 0.001;
				u1 = numerov (V, E1, rmin, rmax, N);
			endwhile
		else
			printf ("Error!\nCounter = %d\n", counter);
			break;
		endif

		counter++;
	endwhile

	E0 = E1;
	u  = u1;

	return;
endfunction
