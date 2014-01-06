function [B_new, L_new] = pushB (B, L, j, b, N)
	Q = defactorize (B, L, j, b, N);

	[U, D, V] = svd (Q, 1);

	[nr, nc] = size (U);

	if j >= 2
		lambda_inv = diag(1./diag(L{j-1}));
	else
		lambda_inv = 1;
	endif

	% every second line forms B^{(j)}_{s_j}
	% -- so we update B^j
	for k = 1:nr/2
		B_new{1,j}(k,:) = U (2*k-1, :);
		B_new{2,j}(k,:) = U (2*k, :);
	end

	B_new{1,j} = lambda_inv * B_new{1,j};
	B_new{2,j} = lambda_inv * B_new{2,j};

	V = ctranspose (V);
	[nr, nc] = size (V);

	% every second column forms B^{(j+1)}_{s_{j+1}}
	% -- so we update B^{j+1}
	for k = 1:nc/2
		B_new{1,j+1}(:,k) = V (:,2*k-1);
		B_new{2,j+1}(:,k) = V (:,2*k);
	end

	% we update the Schmidt coefficients
	% -- in other words, Lambda^j
	L_new{j} = D;

	% all the rest remain the same
	for k = 1:j-1
		B_new{1,k} = B{1,k};
		B_new{2,k} = B{2,k};
		L_new{k} = L{k};
	end

	for k = j+2:2*N
		B_new{1,k} = B{1,k};
		B_new{2,k} = B{2,k};
		L_new{k-1} = L{k-1};
	end

	return;
endfunction
