function pauli = sigma (k)
	matrix = [[0, 1; 1, 0], [0, -i; i, 0], [1, 0; 0,-1]];
	pauli = matrix (:,(2*k-1):2*k);
	return
endfunction
