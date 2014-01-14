function n = altnorm (B, L, N)

	T = matrixT (B, L, N);

	n = T{1};

	for j = 2:2*N
		new_n = n * T{j};
		n = reshape (new_n, size(new_n));
	end

	n = sqrt(norm(n));
	return;

endfunction
