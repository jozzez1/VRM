function v = ground_state (N, periodic_flag)
	H = hamiltonN (N, periodic_flag);
	[V, Lambda] = eigs (H,1,"sa");

	v = V(:,1);
	v = v.';
	return;
endfunction
