function c = mtimes(H, b)
if isfloat(H) && isscalar(H)
      c = HMatrix('pointer', hmatrix_scale (b, H), b.row_cluster, b.col_cluster);
      return;
  end
  
  if isfloat(b) && isscalar(b)
      c = HMatrix('pointer', hmatrix_scale (H, b), H.row_cluster, H.col_cluster);
      return;
  end
 if (isa(b, 'HMatrix'))
    % TODO: We need to check that dimensions and clusters match. 
    Cptr = hmatrix_prod (H, b);
    c = HMatrix('pointer', Cptr, H.row_cluster, b.col_cluster);
    return;
end
  if (ismatrix(b) && ~isscalar(b))
    sz = size(b,2);
    if sz == 1
      c = mvm_hmatrix_avector(H, b);
    else
      c = hmatrix_prod_full(H, b);
    end
    return;
end
  

    
end
