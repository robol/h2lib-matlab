function c = mtimes(H, b)
  if (isa(b, 'HMatrix'))
    % TODO: We need to check that dimensions and clusters match. 
    Cptr = hmatrix_prod (H, b);
    c = HMatrix('pointer', Cptr, H.row_cluster, b.col_cluster);
  else
    sz = size(b);

    % TODO: Handle more columns as well
    if (sz(2) == 1)
      c = mvm_hmatrix_avector(H, b);
    end
  end	   
end
