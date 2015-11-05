function C = uminus(A)
if (isa(A, 'HMatrix'))
      Cptr = hmatrix_uminus(A);
      C = HMatrix('pointer', Cptr, A.row_cluster, A.col_cluster);
  else
      C = -A;
  end
end
