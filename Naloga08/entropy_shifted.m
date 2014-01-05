function E = entropy_shifted (v, N)
% here's the difference -- we fill the matrix in a weird way
P = auxmatrix (N);
for i = 1:2**N
	for j = 1:2**N
		P (i,j) = v(P(i,j));
	end
end

% the rest is the same
[U, D, V] = svd (P);

d = diag (D);

E = 0;
for i = 1:2**N
	ell = d(i)**2;
	E -= ell * log (ell);
end

endfunction
