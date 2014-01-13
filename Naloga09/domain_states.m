function [f1, f2] = domain_states (N, complex_flag)
	f1 = rand_state (N, complex_flag);

	k = 2**N * (-1 + 2**N)+1;
	norma = norm(f1(k));

	f1 = f1 ./norma;
	f2 = f1;
	f1(k) = 2*f1(k);
	
	%test
	f1 - f2

	return;
endfunction
