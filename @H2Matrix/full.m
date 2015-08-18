function H = full(P)
  sz = matrix_size(P)
  n = sz(2);
  m = sz(1);

  H = zeros(sz);

  for i = 1 : n
    H(:,i) = P * [ zeros(i-1,1) ; 1 ; zeros(n-i,1) ];
    % H(:,i) = mvm_hmatrix_avector(P, [ zeros(i-1,1) ; 1 ; zeros(n-i,1) ]);
  end
end
