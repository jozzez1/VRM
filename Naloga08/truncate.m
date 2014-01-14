function [Bnew, Lnew] = truncate (B, L, N, M)
	Bnew = B;
	Lnew = L;

	if M > 1
		for j = 1:2*N-1
			[nr, nc] = size(L{j});

			if nr > M
				Lnew{j} = diag(diag(L{j})(1:M));
				Bnew{1,j} = B{1,j}(:,1:M);
				Bnew{2,j} = B{2,j}(:,1:M);
				Bnew{1,j+1} = B{1,j+1}(1:M,:);
				Bnew{2,j+1} = B{1,j+1}(1:M,:);
			endif
		end
	endif

	return;
endfunction
