function H = hmatrix_full(P)
  sz = hmatrix_size(P);
  n = sz(2);
  m = sz(1);

  H = zeros(sz);

  for i = 1 : n
    H(:,i) = mvm_hmatrix_avector(P, [ zeros(i-1,1) ; 1 ; zeros(n-i,1) ]);
  end
end
