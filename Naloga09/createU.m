function U = createU (b)
	U = [exp((-2)*b), 0, 0, 0;
		0, cosh(2*b), (-1)*sinh(2*b), 0;
		0, (-1)*sinh(2*b), cosh(2*b), 0;
		0, 0, 0, exp((-2)*b)];

	U *= exp(b);

	return;
endfunction
