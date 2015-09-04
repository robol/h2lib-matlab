function c = mtimes(H, b)
  if (isa(b, 'H2Matrix'))
    % TODO: We need to check that dimensions and clusters match. 
    Cptr = h2matrix_prod (H, b);
    c = H2Matrix('pointer', Cptr, H.row_cluster, b.col_cluster);
  else
    sz = size(b);

    % TODO: Handle more columns as well
    if (sz(2) == 1)
      c = mvm_h2matrix_avector(H, b);
    end
  end	   
end
