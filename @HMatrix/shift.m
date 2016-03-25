function B = shift (H, r)
%SHIFT Shift the matrix H and set B = H + rI. 

  B = HMatrix('pointer', hmatrix_shift(H, r), H.row_cluster, H.col_cluster);

end
