function T = matrixT (B, L, N)

	T{1} = kron(conj(B{1,1}), B{1,1}) + kron(conj(B{2,1}), B{2,1});

	for j = 2:2*N
		A1 = L{j-1} * B{1,j};
		A2 = L{j-1} * B{2,j};
		T{j} = kron(conj(A1), A1) + kron(conj(A2), A2);
	end

endfunction
