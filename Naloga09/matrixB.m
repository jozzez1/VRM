function K = matrixB (B, L, k, j, b, N)
	% we are looking at k-th component of the vector, and
	% are grouping j-th and (j+1)-st spin together
	% our chain has 2*N-spins

	U = matrixU (b);

	n (1,1) = dec2bin(k-1, 2*N)(1,j);
	n (1,2) = dec2bin(k-1, 2*N)(1, mod(j,2*N) + 1);
	x = U (bin2dec(n)+1,:);

	if j <= 2*N - 1
		lambda = L{j};
	elseif j == 2*N
		lambda = 1;
	endif

	K = 0;
	for i = 1:4
		s1 = bin2dec(dec2bin (i-1,2)(1,1)) + 1;
		s2 = bin2dec(dec2bin (i-1,2)(1,2)) + 1;

		K += x(i) * B{s1, j} * lambda * B{s2,mod(j, 2*N) + 1};
	end

	return;
endfunction
