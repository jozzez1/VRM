function S = altentropy (L)
	lambda = diag (L);
	[nr, nl] = size(L);

	S = 0;
	for i = 1:nr
		ell = lambda(i)^2;
		S -= ell * log(ell);	
	end
	return;
endfunction
