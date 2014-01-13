function [Bnew, Lnew] = evenstep (U, B, L, N)
	Bnew = B;
	Lnew = L;
	for j = 1:N-1
		[Bnew, Lnew] = pushB (U, 2*j, Bnew, Lnew);
	end
	
	return;
endfunction
