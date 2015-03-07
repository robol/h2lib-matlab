function M = inv(H)
    M = HMatrix('pointer', hmatrix_inv(H));
end