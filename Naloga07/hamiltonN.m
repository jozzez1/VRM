function H = hamiltonN (N, periodic_flag)
	% periodic flag should be either 0 or 1!
	
	% two-particle hamiltonian we obtain from Pauli matrices
	h = zeros (4);
	for i = 1:3
		h = h + kron(sigma(i), sigma(i)); 
	end

	% periodic boundary correction	
	H = zeros (4**N);
	for i = 1:3
		H = H + periodic_flag * kron (sigma(i), eye(4**(N-1)), sigma(i));
	end

	% and part for open boundaries
	for i = 1:(2*N-1)
		H = H + kron (eye(2**(i-1)), h, eye(2**(2*N-1-i)));
	end

	return;
endfunction
