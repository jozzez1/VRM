function [Bnew, Lnew] = fullstepS4 (Uhalf1, Ufull1, Uhalf0, Ufull0, B, L, N)
	Bnew = B;
	Lnew = L;

	[Bnew, Lnew] = fullstep (Uhalf1, Ufull1, Bnew, Lnew, N);
	[Bnew, Lnew] = fullstep (Uhalf0, Ufull0, Bnew, Lnew, N);
	[Bnew, Lnew] = fullstep (Uhalf1, Ufull1, Bnew, Lnew, N);

	return;
endfunction
