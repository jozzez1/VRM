function [Bnew, Lnew] = part_step (B, L, b, N, odd_even_flag)

	for j = 1:N
		if strcmp (odd_even_flag, 'odd')
			k = 2*j - 1;
		elseif strcmp (odd_even_flag, 'even')
			k = 2*j;
		endif

		[Bnew,Lnew] = pushB (B, L, k, b, N);
	end

	return;
endfunction
