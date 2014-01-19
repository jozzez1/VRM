function [E0, u] = eigsground (V, rmax, N)
	H = (-0.5)*mLaplace (rmax, N) + diag(V);

	[u, E0] = eigs (H, 4, "sa");

	return;
endfunction
