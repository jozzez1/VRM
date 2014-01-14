function [E0, d] = regression (logN, b)
	% input logN and b are both row-vectors
	[nr, nc] = size(b);
	X = eye(nc, 2);
	X(:,1) = (-1) .* b.';

	% logN must be a column-vector for this to work
	logN = logN.';

	v = inv(transpose(X) * X) * transpose(X)*logN;

	E0 = v(1);
	d = v(2);

	return;
endfunction
