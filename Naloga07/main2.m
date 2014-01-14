% matlab/octave script for the MPA
% WARNING -- Postponed indefenitely

% initialize the outpu files
fout1 = fopen ("2-periodic.txt", "w+");
fout2 = fopen ("2-nonperiod.txt", "w+");
fout3 = fopen ("2-gauss-rand.txt", "w+");

% calculate the data
printf ("Calculating ... ");
for k = 2:6
	% first we create our eigenvectors
	vo = ground_state (k, 0);
	vp = ground_state (k, 1);
	vr = rand_state (k);

	% now comes the MPA
	[Ao, Lo] = matricize (vo);
	[Ap, Lp] = matricize (vp);
	[Ar, Lr] = matricize (vr);

	% and now we calculate the entropy
	for a = 1:2*k-1
		Ep (k-1,a) = altentropy (Lp{a});
		Eo (k-1,a) = altentropy (Lo{a});
	end

	[alpha(k-1), normdiff(k-1)] = check (Ar, vr);
end
printf ("Done!\n");

% now we dump this to some files
printf ("Dumping results ... ");
for k = 2:6
	fprintf (fout1, "%f %f ", 0.5/k, Ep(k,1));
	fprintf (fout2, "%f %f ", 0.5/k, Eo(k,1));

	% 1st column are open boundary conditions -- both non-compact
	% 2nd column is for periodic boundaries
	% 3rd column is for compact case -- open
	% 4th column is for compact case -- periodic
	fprintf (fout3, "%d %f %f %f %f\n", 2*k, Eno(k), Enp(k), Eo(k,k), Ep(k,k));
end
