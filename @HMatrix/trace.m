function t=trace(H)
if (isa(H, 'HMatrix'))
      t=hmatrix_trace(H);
end
end
