function [Smax, S] = gaus_entropy (N)
	v = rand_state (N);
	
	[A, L] = matricize (v);
	
	fout = fopen ("nakljuc.txt", "w+");

	for i = 1:2*N-1
		a (i) = 0.5 * i/N;
		S (i) = altentropy (L{i});
		Smax (i) = maxentropy (N, i);

		fprintf (fout, "%f\t%f\t%f\n", a(i), S(i), Smax(i));
	end

	fclose (fout);
	return;
endfunction
