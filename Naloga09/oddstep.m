function [Bnew, Lnew] = oddstep (U, B, L, N, M)
	Bnew = B;
	Lnew = L;
	for j=1:N
		[Bnew, Lnew] = pushB (U, 2*j-1, Bnew, Lnew, M);
	end

	return;
endfunction
