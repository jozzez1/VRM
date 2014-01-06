function Q = defactorize (B, L, j, b, N)
	
	for k = 1:4
		if j >= 2
			K{k} = L{j-1} * matrixB (B, L, k, j, b, N);
		else
			K{k} = matrixB (B, L, k, 1, b, N);
		endif
	end

	[nr, nc] = size (K{1});
	Q = zeros (2*nr, 2*nc);

	Q(1:nr, 1:nc) = K{1};
	Q(1:nr, nc+1 : 2*nc) = K{2};
	Q(nr+1:2*nr, 1:nc) = K{3};
	Q(nr+1:2*nr, nc+1:2*nc) = K{4};

	return;
endfunction
