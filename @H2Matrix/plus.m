function C = plus(A,B)
  if (isa(B, 'H2Matrix'))
      Cptr = h2matrix_sum (A, B);
      C = H2Matrix('pointer', Cptr, A.row_cluster, A.col_cluster);
  else
      C = full(A) + B;
  end
end
