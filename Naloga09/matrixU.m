function U = matrixU (b)
	U = [exp(-2*b), 0, 0, 0; 0, cosh(2*b), -sinh(2*b), 0;
		0, -sinh(2*b), cosh(2*b), 0; 0, 0, 0, exp(-2*b)];

	return;
endfunction

