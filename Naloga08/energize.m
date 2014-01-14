function [E0, B, L] = energize (v, bmax, number_of_steps, M)
	[nr, nc] = size(v);
	N = log2(nc)/2;

	% M is truncation parameter, M > 1

	[B, L] = matricize (v);

	bstep = bmax/number_of_steps;

	% we will use complex method S3 -- Trotter-Suzuki
	% we do that, because we have dissipation
	p1 = (1 + i/cbrt(3))/8;
	p2 = 2*p1;
	p3 = 1/4;
	p4 = conj(p2);
	p5 = conj(p1);
	
	% matrices U
	Uhalf_p1 = createU (p1 * bstep/2);
	Ufull_p1 = createU (p1 * bstep);
	Uhalf_p2 = createU (p2 * bstep/2);
	Ufull_p2 = createU (p2 * bstep);
	Uhalf_p3 = createU (p3 * bstep/2);
	Ufull_p3 = createU (p3 * bstep);
	Uhalf_p4 = createU (p4 * bstep/2);
	Ufull_p4 = createU (p4 * bstep);
	Uhalf_p5 = createU (p5 * bstep/2);
	Ufull_p5 = createU (p5 * bstep);

	report = fopen ("report.dat", "w+");
	report1= fopen("report1.dat", "w+");

	K = ceil(5/bstep);

	for k = 1:number_of_steps

		[B, L] = fullstep (Uhalf_p1, Ufull_p1, B, L, N, M);
		[B, L] = fullstep (Uhalf_p2, Ufull_p2, B, L, N, M);
		[B, L] = fullstep (Uhalf_p3, Ufull_p3, B, L, N, M);
		[B, L] = fullstep (Uhalf_p4, Ufull_p4, B, L, N, M);
		[B, L] = fullstep (Uhalf_p5, Ufull_p5, B, L, N, M);

		norma = altnorm (B, L, N);
		bet = bstep * k;
		fprintf (report, "%f\t%f\n", bet, log(norma));
		if bet >= 5;
			b(k+1-K) = bet;
			logN(k+1-K) = log(norma);
			fprintf (report1, "%f\t%f\n", b(k+1-K), logN(k+1-K));
		endif
	end
	fclose (report);
	fclose (report1);

	[E0, d] = regression (logN, b);

	return;
	
endfunction
