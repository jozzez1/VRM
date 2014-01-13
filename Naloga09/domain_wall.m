function domain_wall (N, tmax, steps, M)

	[v1, v2] = domain_states (N, 1);

	[B1, L1] = matricize (v1, M);
	[B2, L2] = matricize (v2, M);

	domain = fopen ("domain.txt", "w+");
	prof = fopen ("profile.txt", "w+");

	tstep = tmax/steps;
	Uf_half = createU ((-i)*tstep/2);
	Uf_full = createU ((-i)*tstep);

	for k = 0:steps
		s11 = sigma_j (B1, L1, B1, L1, N);
		s22 = sigma_j (B2, L2, B2, L2, N);
		s12 = sigma_j (B1, L1, B2, L2, N);

		fprintf (domain, "%f ", k*tstep);
		for j = 1:2*N
			fprintf (domain, "%f ", s11(j) + s22(j) - 2*real(s12(j)));
			fprintf (prof, "%f\t%f\n", k*tstep, s11(j) + s22(j) - 2*real(s12(j)));
		end
		fprintf (domain, "\n");
		fprintf (prof, "\n");

		[B1, L1] = fullstep (Uf_half, Uf_full, B1, L1, N, M);
		[B2, L2] = fullstep (Uf_half, Uf_full, B2, L2, N, M);
	end
	fclose (domain);
endfunction
