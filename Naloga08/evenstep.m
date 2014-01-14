function [Bnew, Lnew] = evenstep (U, B, L, N, M)
	Bnew = B;
	Lnew = L;
	for j = 1:N-1
		[Bnew, Lnew] = pushB (U, 2*j, Bnew, Lnew, M);
	end
	
	return;
endfunction
