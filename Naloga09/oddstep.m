function [Bnew, Lnew] = oddstep (U, B, L, N)
	Bnew = B;
	Lnew = L;
	for j=1:N
		[Bnew, Lnew] = pushB (U, 2*j-1, Bnew, Lnew);
	end

	return;
endfunction
