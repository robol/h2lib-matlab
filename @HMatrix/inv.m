function M = inv(H)
    M = HMatrix('pointer', hmatrix_inv(H), H.row_cluster, H.col_cluster);
end