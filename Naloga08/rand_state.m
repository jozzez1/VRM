function v = rand_state (N)
	% we create two row-vectors
	v_real = normrnd (0, 3, 1, 4**N);
	v_imag = normrnd (0, 3, 1, 4**N);

	% combine them into one complex vector
	v = v_real + i*v_imag;

	% let's normalize it to filter away the
	% size effects
	v = v ./ norm (v);

	return;
endfunction
