function M = mldivide(H, A)
    if (isa(A, 'HMatrix'))
        Hinv = HMatrix('pointer', hmatrix_inv (H));
        M = Hinv * A;
    else
        % TODO: Implement this as a linear system solver. 
        M = full(H) \ A;
    end
end
