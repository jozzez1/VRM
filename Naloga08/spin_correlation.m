function C = spin_correlation (B, L, N, M)
	% B, L correspond to a unnormalized state -- we have to normalize it
	v = reconstruct (B, L, N);
	v = v ./ norm(v);
	
	% now we take the new matrices B, L
	[B, L] = matricize (v);

	% and now's the time for the real fun
	T = matrixT (B, L, N);

	% we create the spin matrices
	for j = 2:2*N-1
		V{j} = kron(conj(L{j-1}*B{1,j}), L{j-1}*B{1,j}) - kron(conj(L{j-1}*B{2,j}), L{j-1}*B{2,j});
	end

	for j = 2:2*N-1
		for k = 2:2*N-1
			Cjk = T{1};
			for i = 2:2*N
				if i == j
					mat = V{j};
				elseif i == k
					mat = V{k};
				else
					mat = T{i};
				endif
				temp = Cjk * mat;
				Cjk = reshape(temp, size(temp));
			end
			Ca(j,k) = Cjk;
		end
	end

	C = Ca (2:2*N-1,2:2*N-1);

	return;
endfunction
