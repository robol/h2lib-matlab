function c = mtimes(H, b)
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
