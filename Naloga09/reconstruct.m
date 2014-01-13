function v = reconstruct (B, L, N)	
	for i = 1:4**N
		n = dec2bin (i-1, 2*N);
		T = 1;
		for j = 1:2*N
			sj = bin2dec(n(1,j)) + 1;

			lambda = 1;
			if j > 1
				lambda = L{j-1};
			end

			b = T * lambda * B{sj, j};
			T = reshape (b, size(b));
		end
		v (i) = T;
		clear b T;
	end

	return;
endfunction
