function s = sigma_j (B1, L1, B2, L2, N)

	for j = 1:2*N
		if j == 1
			Lam1 = 1;
			Lam2 = 1;
		else
			Lam1 = L1{j-1};
			Lam2 = L2{j-1};
		endif
		T{j} = kron(conj(Lam1*B1{1,j}), Lam2*B2{1,j}) + kron(conj(Lam1*B1{2,j}), Lam2*B2{2,j});
		V{j} = kron(conj(Lam1*B1{1,j}), Lam2*B2{1,j}) - kron(conj(Lam1*B1{2,j}), Lam2*B2{2,j});
	end

	for j = 1:2*N
		if j == 1
			sj = V{1};
		else
			sj = T{1};
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
