function c = mtimes(H, b)

  if isfloat(H) && isscalar(H)
      n = matrix_size(b);
      n = n(1);
      D = HMatrix('tridiagonal', b.row_cluster, b.col_cluster, ...
                  H * ones(n, 1), zeros(n-1,1), zeros(n-1,1));
      c = mtimes(D, b);
      return;
  end
  
  if isfloat(b) && isscalar(b)
      n = matrix_size(H);
      n = n(1);
      D = HMatrix('tridiagonal', H.row_cluster, H.col_cluster, ...
                  b * ones(n, 1), zeros(n-1,1), zeros(n-1,1));
      c = mtimes(D, H);
      return;
  end

  if (isa(b, 'HMatrix'))
    % TODO: We need to check that dimensions and clusters match. 
    Cptr = hmatrix_prod (H, b);
    c = HMatrix('pointer', Cptr, H.row_cluster, b.col_cluster);
  else
    
    sz = size(b,2);
    c=zeros(size(b,1),sz);
    % TODO: Handle more columns as well
    for j=1:sz
      c(:,j) = mvm_hmatrix_avector(H, b(:,j));
    end
  end	   
end
