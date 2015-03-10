ps = 1;

% try 

addpath('..');

fprintf (ps, 'Testing the main operations of H2Lib\n\n')

if ~exist('n', 'var')
    n = 100; 
end

a = randn(n,1); b = randn(n-1,1); c = randn(n-1,1);
a2 = randn(n,1); b2 = randn(n-1,1); c2 = randn(n-1,1);
cc = Cluster(n, 4);

H = inv(HMatrix('tridiagonal', cc, cc, a, b, c));
H2 = HMatrix('tridiagonal', cc, cc, a2, b2, c2);

fH = full(H);
fH2 = full(H2);

tic; 
H3 = H * H2;
r = toc;

tic;
fH3 = fH * fH2;
rf = toc;

res = norm(full(H3) - fH3);
fprintf (ps, ' > Residue of full(H*H2) - full(H)*full(H2) = %e\n', res);
fprintf (ps, ' > Time needed for structured op: %e\n', r);
fprintf (ps, ' > Time needed for full op: %e\n\n', rf);

tic;
H3 = H + H2;
r = toc;

tic;
fH3 = fH + fH2;
rf = toc;

res = norm(fH3 - full(H3));

fprintf (ps, ' > Residue of full(H + H2) - full(H) - full(H2) = %e\n', res);
fprintf (ps, ' > Time needed for structured op: %e\n', r);
fprintf (ps, ' > Time needed for full op: %e\n\n', rf);

tic;
H3 = H \ H2;
r = toc;

tic;
fH3 = fH \ fH2;
rf = toc;

res = norm(fH3 - full(H3));
fprintf (ps, ' > Residue of full(H + H2) - full(H) \\ full(H2) = %e\n', res);
fprintf (ps, ' > Time needed for structured op: %e\n', r);
fprintf (ps, ' > Time needed for full op: %e\n\n', rf);

tic;
H3 = H / H2;
r = toc;

tic;
fH2 = fH / fH2;
rf = toc;
res = norm(fH3 - full(H3));
fprintf (ps, ' > Residue of full(H + H2) - full(H) / full(H2) = %e\n', res);
fprintf (ps, ' > Time needed for structured op: %e\n', r);
fprintf (ps, ' > Time needed for full op: %e\n\n', rf);

tic;
H3 = inv(H);
r = toc;

tic;
fH3 = inv(fH);
rf = toc;

res = norm(full(H3) * fH - eye(n));
fprintf (ps, ' > Residue of full(H * inv(H)) - I = %e\n', res);
fprintf (ps, ' > Time needed for structured op: %e\n', r);
fprintf (ps, ' > Time needed for full op: %e\n\n', rf);

fprintf (ps, '\n');

% catch e
%     fprintf (ps, 'Problem running the tests, please check them manully\n');
% end;
