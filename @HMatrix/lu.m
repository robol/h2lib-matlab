function [L,U]=lu(H)
[L,U]=hmatrix_LU(H);
L= HMatrix('pointer', L, H.row_cluster, H.col_cluster);
U= HMatrix('pointer', U, H.row_cluster, H.col_cluster);
end
