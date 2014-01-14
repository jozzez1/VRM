function [ s ]=schmidt( vf,n,k )

  matrika=sparse(2^k, 2^(n-k), 0);
  
  %sedaj mesto v seznamu pretvorimo v binarno kodo in vzamemo prvih k in drugih n-k
  %oboje pretvorimo nazaj v integrer in dobimo vrstico, stolpec :-)
   
  for i=1:2^n
    
    j=i-1;
    
    j_bin=dec2bin(j, n);
    
    k_part_bin=j_bin(1:k);
    rest_bin=j_bin((k+1):n);
    
    vrstica=bin2dec(k_part_bin)+1;
    stolpec=bin2dec(rest_bin)+1;
     
    matrika(vrstica, stolpec)=vf(i);
    
  end
  
  %sedaj bomo nardili SVD razcep
  
  s=svds(matrika, 1000);
  
endfunction
