function [Bnew, Lnew] = fullstep (Uhalf, Ufull, B, L, N)
	Bnew = B;
	Lnew = L;

	[Bnew, Lnew] = oddstep (Uhalf, Bnew, Lnew, N);
	[Bnew, Lnew] = evenstep (Ufull, Bnew, Lnew, N);
	[Bnew, Lnew] = oddstep (Uhalf, Bnew, Lnew, N);

	return;
endfunction
