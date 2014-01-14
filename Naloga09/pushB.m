function [Bnew, Lnew] = pushB (mU, j, B, L, M)
	% first we create the matrices Bjj (for interaction)
	for k = 1:4
		Bjj {k} = createBnew (mU, k, j, B, L);
	end

	if j == 1
		lambda = 1;
	elseif j >= 2
		lambda = diag(L{j-1});
	endif

	% now we create the matrix Q
	[nr, nc] = size(Bjj{1});
	Q = zeros (2*nr, 2*nc);

	for i = 1:nr
		for k = 1:nc
			Q(i,k) = lambda(i) * Bjj{1}(i,k);
			Q(i,k+nc) = lambda(i) * Bjj{2}(i,k);
			Q(i+nr, k) = lambda(i) * Bjj{3}(i,k);
			Q(i+nr, k+nc) = lambda(i) * Bjj{4}(i,k);
		end
	end

	[U,D,V] = svd (Q,1);
	V = ctranspose(V);

	% we grab the new matrices B and L
	Bnew = B;
	Lnew = L;

%	Linv = diag(1 ./ lambda);
	Linv = inv(diag(lambda));

	[nr,nc] = size(U);
	Bnew{1,j} = Linv * U(1:nr/2,:);
	Bnew{2,j} = Linv * U(nr/2 + 1:nr,:);

	[nr,nc] = size(V);
	Bnew{1,j+1} = V(:,1:nc/2);
	Bnew{2,j+1} = V(:,nc/2 + 1:nc);

	Lnew{j} = D;

	% and here we truncate
	[nr, nc] = size(Lnew{j});
	if ((nc > M) && (M > 1))
		Lnew{j} = diag(diag(Lnew{j})(1:M));
		Bnew{1,j} = Bnew{1,j}(:,1:M);
		Bnew{2,j} = Bnew{2,j}(:,1:M);
		Bnew{1,j+1} = Bnew{1,j+1}(1:M,:);
		Bnew{2,j+1} = Bnew{2,j+1}(1:M,:);
	endif

	return;
endfunction
