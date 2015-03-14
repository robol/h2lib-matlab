function C = plus(A,B)
  if (isa(B, 'HMatrix'))
      Cptr = hmatrix_sum (A, B);
      C = HMatrix('pointer', Cptr, A.row_cluster, A.col_cluster);
  else
      C = full(A) + B;
  end
end
