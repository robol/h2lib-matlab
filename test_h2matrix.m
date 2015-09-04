n = 16;

a = randn(n,1);
b = randn(n,1);
c = randn(n,1);
% cc = Cluster(n, 4);
cc = Cluster(n, 4);
H = H2Matrix('tridiagonal', cc, cc, a, b, c);
full(H)

H2 = H2Matrix('tridiagonal', cc, cc, a + 1, b, c);
fprintf ('Residue of sum: %e\n', norm (full(H + H2) - full(H) - full(H2)));
fprintf ('Residue of prod: %e\n', norm (full(H * H2) - full(H) * full(H2)));
fprintf ('Residue of mldivide: %e\n', norm (full(H \ H2) - full(H) \ full(H2)));
