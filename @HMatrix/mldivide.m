function M = mldivide(H, A)
    if (isa(A, 'HMatrix'))
        l = hmatrix_mldivide(H,A);
        M = [];
        for i = 1 : size(A,2)
          M = [ M, HMatrix('pointer', l(i)) ];
        end
    else
        % TODO: Implement this as a linear system solver. 
        M = full(H) \ A;
    end
end
