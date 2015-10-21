n = 1280;

if ~exist('avoid_redefinition', 'var') || avoid_redefinition == false

a = randn(n,1);
b = randn(n,1);
c = randn(n,1);
% cc = Cluster(n, 4);
cc = Cluster(n, 96);
H = H2Matrix('tridiagonal', cc, cc, a, b, c);
H2 = H2Matrix('tridiagonal', cc, cc, a + 1, b, c);

end

fprintf ('Residue of sum: %e\n', norm (full(H + H2) - full(H) - full(H2)));
fprintf ('Residue of prod: %e\n', norm (full(H * H2) - full(H) * full(H2)));
fprintf ('Residue of mldivide: %e\n', norm (full(H \ H2) - full(H) \ full(H2)));

return

% HF = full(H);

[LL,UU] = lu(sparse(full(H)), 0.0); LL = full(LL); UU = full(UU);
%l = H \ [ H2 H2 ]; L = l(1); R = l(2);
%fprintf ('Residue: %e\n', norm(full(L) - LL));
%fprintf ('Residue: %e\n', norm(full(R) - UU));

l = H \ [ H2 H2 ]; L = l(1); U = l(2);
fprintf ('Residue: %e\n', norm(full(L) - UU \ full(H2)));
% fprintf ('Residue: %e\n', norm(full(L) - (UU \ full(H2))));
% fprintf ('Residue: %e\n', norm(full(R) - (UU \ full(H2))));


% fprintf ('Has H changed?: %e\n', norm(full(H) - HF));
