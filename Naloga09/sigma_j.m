function s = sigma_j (B, L, N)

	T = matrixT (B, L, N);

	for j = 1:2*N
		if j == 1
			Lam = 1;
		else
			Lam = L{j-1}
		endif

		V{j} = kron(conj(Lam * B{1,j}), Lam * B{1,j}) - kron(conj(Lam*B{2,j}), Lam*B{2,j});
	end

	for j = 1:2*N
		if j == 1
			sj = V{1};
		else
			sj = T{1}
		endif

		for k = 2:2*N
			if k == j
				mat = V{j};
			else
				mat = T{k};
			endif

			temp = sj * mat;
			sj = reshape (temp, size(temp));
		end
		s(j) = sj;
	end

	return;
endfunction
