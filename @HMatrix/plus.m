function C = plus(A,B)
  if (isa(B, 'HMatrix'))
      Cptr = hmatrix_sum (A, B);
      C = HMatrix('pointer', Cptr);
  else
      C = full(A) + B;
  end
end
