function [alpha, diffnorm] = check (A,v)
	% with this function we check wheather or not we get the
	% same result:
	% ==============
	% alpha    -- angle between two vectors
	% normdiff -- difference in norms of the two
	%
	% the reconstructed vector from matrices 'A' will be called 'x'

	[nr, nc] = size (v);
	N = log2 (nc)/2;
	for i = 1:4**N
		T	= 1;
		k	= dec2bin (i-1, 2*N);
		for j = 1:2*N
			b = T * A{bin2dec(k(1,j))+1, j};
			T = reshape(b, size(b));
		end
		x (i) = T;

		clear b T;
	end

	alpha	= norm(acos ((x * v')/(norm(x) * norm(v))));
	diffnorm= norm (x - v);

	return;
endfunction
