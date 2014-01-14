function B = basis2x2 (i)
	A = [[1,0;0,0],[0,1;0,0],[0,0;1,0],[0,0;0,1]];
	B = A (:,2*i-1:2*i);
	return;
endfunction
