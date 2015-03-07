function c = mtimes(H, b)
  sz = size(b);
  
  if (sz(2) == 1)
      % Matrix vector product
      c = mvm_hmatrix_avector(H, b);
  else
      % Matrix matrix product. We consider the only case where b
      % is a n x n matrix with the same dimensions of H, as of now. 
      if (isa(b, 'HMatrix'))
          % TODO: We need to check that dimensions and clusters match. 
          Cptr = hmatrix_prod (H, b);
          c = HMatrix('pointer', Cptr);
      end
  end
end
