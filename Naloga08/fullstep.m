function [Bnew, Lnew] = fullstep (Uhalf, Ufull, B, L, N, M)
	Bnew = B;
	Lnew = L;

	[Bnew, Lnew] = oddstep (Uhalf, Bnew, Lnew, N, M);
	[Bnew, Lnew] = evenstep (Ufull, Bnew, Lnew, N, M);
	[Bnew, Lnew] = oddstep (Uhalf, Bnew, Lnew, N, M);

	return;
endfunction
