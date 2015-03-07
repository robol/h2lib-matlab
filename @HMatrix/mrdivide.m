function M = mrdivide(H, A)
    if (isa(A, 'HMatrix'))
        Ainv = HMatrix('pointer', hmatrix_inv(A));
        M = H * Ainv;
        release(Ainv);
    else
        % TODO: Implement this as a linear system solver
        M = A / full(H);
    end
end