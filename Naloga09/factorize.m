function [A, L, P] = factorize (V)
	% we start by rearraging the components of 'V' into a new matrix 'X'

	[nr, nc] = size(V);

	for i = 1:nr
		for j = 1:nc
			X (2*i - 1 + floor(j/(1 + nc/2)), 1 + mod(j-1, nc/2)) = V(i,j);
		end
	end

	% now we make an SVD
	[A,L,P] = svd (X,1);

	return;
endfunction

