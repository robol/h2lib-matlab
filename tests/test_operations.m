ps = 2;

try; 

addpath('..');

fprintf (ps, 'Testing the main operations of H2Lib\n\n')

n = 100; 
a = randn(n,1); b = randn(n-1,1); c = randn(n-1,1);
a2 = randn(n,1); b2 = randn(n-1,1); c2 = randn(n-1,1);
cc = Cluster(n, 4);

H = inv(HMatrix('tridiagonal', cc, cc, a, b, c));
H2 = HMatrix('tridiagonal', cc, cc, a2, b2, c2);

H3 = H * H2;

res = norm(full(H3) - full(H) * full(H2));
fprintf (ps, ' > Residue of full(H*H2) - full(H)*full(H2) = %e\n', res);

H3 = H + H2;
res = norm(full(H3) - full(H) - full(H2));
fprintf (ps, ' > Residue of full(H + H2) - full(H) - full(H2) = %e\n', res);

H3 = H \ H2;
res = norm(full(H3) - full(H) \ full(H2));
fprintf (ps, ' > Residue of full(H + H2) - full(H) \\ full(H2) = %e\n', res);

H3 = H / H2;
res = norm(full(H3) - full(H) / full(H2));
fprintf (ps, ' > Residue of full(H + H2) - full(H) / full(H2) = %e\n', res);

H3 = inv(H);
res = norm(full(H3) * full(H) - eye(n));
fprintf (ps, ' > Residue of full(H * inv(H)) - I = %e\n', res);

fprintf (ps, '\n');

catch;
end;
