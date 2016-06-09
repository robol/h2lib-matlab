function Hl = leading(H)    
    Hl = HMatrix();
    Hl.row_cluster = Cluster();
    Hl.row_cluster.parent = H.row_cluster;
    Hl.col_cluster = Cluster();   
    Hl.col_cluster.parent = H.col_cluster;
    hmatrix_leading(H, Hl, Hl.row_cluster, Hl.col_cluster);
end