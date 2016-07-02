function Hl = trailing(H,rc,cc)    
    Hl = HMatrix();
    Hl.row_cluster = rc;
    Hl.col_cluster = cc;    
    hmatrix_trailing(H, Hl, Hl.row_cluster, Hl.col_cluster);
end
