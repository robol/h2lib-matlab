function M = mldivide(H, A)
    if (isa(A, 'HMatrix'))
        l = hmatrix_mldivide(H,A);
        M = [];
        for i = 1 : size(A,2)
          M = [ M, HMatrix('pointer', l(i), H.row_cluster, A.col_cluster) ];
        end
    else
	% Efficient H Matrix linear system solver. 
	M = hmatrix_mldivide_full(H, A);
    end
end
