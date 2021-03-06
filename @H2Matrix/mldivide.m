function M = mldivide(H, A)
    if (isa(A, 'H2Matrix'))        
        l = h2matrix_mldivide(H,A);
        M = [];
        for i = 1 : size(A,2)
          M = [ M, H2Matrix('pointer', l(i), H.row_cluster, H.col_cluster) ];
        end
    else
        % TODO: Implement this as a linear system solver. 
        M = full(H) \ A;
    end
end
