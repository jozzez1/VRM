function A = auxmatrix (N)
	% first we create auxilliary matrix
	A = [1,2;3,4];

	for j = 1:(N-1)
		B = zeros (2**(1+j));
		for i = 1:4
			B += kron (basis2x2(i),A .+ (i-1)*4**j);
		end
		A = B;
	end
	return;
endfunction

