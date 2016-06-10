function Hl = trailing(H)    
    Hl = HMatrix();
    Hl.row_cluster = Cluster();
    Hl.row_cluster.parent = H.row_cluster;
    Hl.col_cluster = Cluster();   
    Hl.col_cluster.parent = H.col_cluster;
    hmatrix_trailing(H, Hl, Hl.row_cluster, Hl.col_cluster);
end
