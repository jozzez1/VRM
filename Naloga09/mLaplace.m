function L = mLaplace (rmin, rmax, N)
	h = (rmax - rmin) / N;
	L = 0;	
	L += diag ((-49/18) * ones(1,N));
	L += diag (1.5 * ones(1,N-1), 1);
	L += diag (1.5 * ones(1,N-1),-1);
	L += diag ((-3.0/20) * ones(1,N-2), 2);
	L += diag ((-3.0/20) * ones(1,N-2),-2);
	L += diag ((1.0/90) * ones (1,N-3), 3);
	L += diag ((1.0/90) * ones (1,N-3),-3);
	L /= h**2;

	return;
endfunction
