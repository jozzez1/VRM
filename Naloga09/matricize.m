function [B, L] = matricize (psi)
	% psi is our quantum spin state
	% 2^n = nc => 2 = log_2 (nc)
	[nr, nc] = size(psi);
	n = log2 (nc);

	% now we create families of matrices
	P{1} = psi;
	for i = 1:n - 1
		[U{i}, L{i}, V{i}] = factorize (P{i});
		P{i+1} = L{i} * ctranspose(V{i});

		[nr, nc] = size(U{i});
		
		for j = 1:nr/2
			% spin up matrix A
			A{1,i}(j,:) = U{i}(2*j-1,:);

			% spin down matrix A
			A{2,i}(j,:) = U{i}(2*j,:);
		end
	end

	% last one is a bit different
	A{1,n} = P{n}(:,1);
	A{2,n} = P{n}(:,2);

	% now we extract matrices B from them
	B{1,1} = A{1,1};
	B{2,1} = A{2,1};

	for i = 2:n
		Linv = diag (1 ./ diag(L{i-1}));
		B{1,i} = Linv * A {1,i};
		B{2,i} = Linv * A {2,i};
	end

	clear A U V;

	return;
endfunction
