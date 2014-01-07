function check (B, Bnew, N)

	for i = 1:2*N
		[nr1, nc1] = size(B{1,i});
		[nr2, nc2] = size(Bnew{1,i});

		printf ("dim1 = %dx%d\tdim2 = %dx%d\n", nr1,nc1, nr2,nc2);
	end

endfunction
