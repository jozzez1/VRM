function [Bnew, Lnew] = energy_step (B, L, b, N)

	% Suzuki-Trotter
	[Bnew, Lnew] = part_step (B, L, 0.5 * b, N, 'odd');
	[Bnew, Lnew] = part_step (Bnew, Lnew, b, N, 'even');
	[Bnew, Lnew] = part_step (Bnew, Lnew, 0.5 * b, N, 'odd');
	
	return;
endfunction
