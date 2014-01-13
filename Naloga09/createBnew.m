function Bnew = createBnew (U, k, j, B, L)
	% k = 1,2,3,4
	% k = 1 ---> s1,s2 = 0,0
	% k = 2 ---> s1,s2 = 0,1
	% ... you get the picture

	Bnew = 0;
	% now we make the multiplications
	for S = 1:4
		s1 = bin2dec(dec2bin(S-1,2)(1,1)) + 1;
		s2 = bin2dec(dec2bin(S-1,2)(1,2)) + 1;

		Bnew += U(k,S) * B{s1, j} * L{j} * B{s2, j+1}; 
	end

	return;
endfunction
