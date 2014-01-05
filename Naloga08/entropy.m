function E = entropy (v, N, A)

%------------------------------------
%let's create the matrix for SVD
%------------------------------------
[ratio, P] = bipartitize (v, A, N);

% and now comes the SVD part already ... yay!
[U, D, V] = svd (P, 1);

% we take the diagonal part for entropy calculation
d = diag (D);

% now to calculate entropy/energy
E = 0;
for i = 1:size(d)
	ell = d(i)**2;
	E -= ell * log (ell);
end

return;
endfunction
