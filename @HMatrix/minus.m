function C = minus(A,B)
  if (isa(B, 'HMatrix'))
      Cptr = hmatrix_minus (A, B);
      C = HMatrix('pointer', Cptr);
  else
      C = full(A) - B;
  end
end
