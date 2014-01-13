function domain_wall (N, tmax, steps)
%	v = rand_state (N, 0);
	v = zeros (1,4**N);		% at the beginning, there was nothing
	v ((2**N)*(2**N - 1)) = 1;	% then I create my state
	[B, L] = matricize (v);		% we MPA this state

	domain = fopen ("domain.txt", "w+");

	tstep = tmax/steps;
	Uf_half = createU ((-i)*tstep/2);
	Uf_full = createU ((-i)*tstep);

	for k = 1:steps
		[B, L] = fullstep (Uf_half, Uf_full, B, L, N);
		s = sigma_j (B, L, N);

		fprintf (domain, "%f ", k*tstep);
		for j = 1:2*N
			fprintf (domain, "%f ", s(j));
		end
		fprintf (domain, "\n");
	end
	fclose (domain);
endfunction
