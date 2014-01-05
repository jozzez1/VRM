function [ratio, P] = bipartitize (v, A, N)
	% total means total number of spins in our chain
	% A means actually number of spins in block A
	% ----------------------------
	% total dimension is 2**(2*N)
	%   -> block A -- 2**A
	%   -> block B -- 2**(2*N - A)
	% ----------------------------

	imax = 2**(2*N - A);
	jmax = 2**A;

	ratio = 4**(A - N);

	for i = 1:imax
		for j = 1:jmax
			P (j, i) = v(i + (j-1)*imax);
		end
	end
	return;
endfunction

