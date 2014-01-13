function v = rand_state (N, complex_flag)
	% we create two row-vectors
	v = normrnd (0, 3, 1, 4**N);

	if complex_flag == 1
		v_imag = normrnd (0, 3, 1, 4**N);
		v += i*v_imag;
	endif

	% let's normalize it to filter away the size effects
	v = v ./ norm (v);

	return;
endfunction
