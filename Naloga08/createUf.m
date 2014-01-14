function Uf = createUf (t)

	Ufp = exp(i * (-t)) * cos(2*t) * eye(4);
	Ufm = i*exp(i * (-t))*sin(2*t) * [1, 0, 0, 0; 0, 0, 1, 0; 0, 1, 0, 0; 0, 0, 0, 1];

	Uf = Ufp + Ufm;

	return;
endfunction
