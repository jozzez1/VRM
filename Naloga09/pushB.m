function [B_new, L_new] = pushB (B, L, j, b, N)
	Q = defactorize (B, L, j, b, N);

	[U, D, V] = svd (Q, 1);

	[nr, nc] = size (U);

	% so we avoid Lambda fiascos
	if j == 1
		lambda_inv = 1/sqrt(sumsq( diag(L{j+1})));
	else
		lambda_inv = diag(1./diag(L{j-1}));
	endif

	% another Lambda fiasco -- we have to preserve the Schmidt rank
	Mj = 2**(N - abs(j - N));

	% every second line forms B^{(j)}_{s_j}
	% -- so we update B^j
	for k = 1:nr/2
		B_new{1,j}(k,:) = U (2*k-1, 1:Mj);	% here we truncate to the Schmidt rank
		B_new{2,j}(k,:) = U (2*k, 1:Mj);	% truncate to the Schmidt rank
	end

	B_new{1,j} = lambda_inv * B_new{1,j};		% truncate to the Schmidt rank
	B_new{2,j} = lambda_inv * B_new{2,j};		% truncate to the Schmidt rank

	V = ctranspose (V);
	[nr, nc] = size (V);

	% every second column forms B^{(j+1)}_{s_{j+1}}
	% -- so we update B^{j+1}
	for k = 1:nc/2
		B_new{1,j+1}(:,k) = V (1:Mj,2*k-1);	% yup, here is truncation too
		B_new{2,j+1}(:,k) = V (1:Mj,2*k);	% truncation
	end

	% we update the Schmidt coefficients
	% -- in other words, Lambda^j
	L_new{j} = diag(diag(D)(1:Mj));			% ok, that's the last of truncation, I promise :D

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
