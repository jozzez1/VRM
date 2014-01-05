% matlab/octave script for homework
%

% initialize the files 
fout1 = fopen ("1-periodic.txt", "w+");		% periodic hamiltonian
fout2 = fopen ("1-nonperiod.txt", "w+");	% aperiodic hamiltonian
fout3 = fopen ("1-non-and-compact.txt", "w+");	% every second spin

% calculate the data
printf ("Calculating ... ");
for k = 2:6
	% first we create our eigenvectors
	vo = rand_state (k, 0);
	vp = rand_state (k, 1);

	% then we calculate entropies for different
	% bi-partitions
	for a = 1:2*k-1
		Ep (k,a) = entropy (vp, k, a);
		Eo (k,a) = entropy (vo, k, a);
	end
	
	% and for different groupings
	Enp (k) = entropy_shifted (vp, k);
	Eno (k) = entropy_shifted (vo, k);
end
printf ("Done!\n");

% make the output to the respective files
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

fprintf (fout1, "\n");
fprintf (fout2, "\n");

for a = 2:11
	for k = 2:6
		fprintf (fout1, "%f %f ", 0.5*a/k, Ep(k,a));
		fprintf (fout2, "%f %f ", 0.5*a/k, Eo(k,a));
	end
	fprintf (fout1, "\n");
	fprintf (fout2, "\n");
end

fclose (fout1);
fclose (fout2);
fclose (fout3);

printf ("Done!\n");
