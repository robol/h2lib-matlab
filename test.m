%A11 = rand(4); 
%A22 = rand(4,4); 
%a = rand(4,1); 
%b = rand(4,1); 
%c = rand(4,1); 
%d = rand(4,1);

%A = hmatrix_create(A11,A22,a,b,c,d);
%x = mvm_hmatrix_avector(A, [ zeros(7,1); 1 ]);

n = 10;
a = randn(n,1);
b = randn(n-1,1);
c = randn(n-1,1);

d = randn(n,1);
ee = randn(n-1,1);
f = randn(n-1,1);

T1 = diag(a) + diag(b,1) + diag(c,-1);
T2 =  diag(d) + diag(ee,1) + diag(f,-1);

P1 = hmatrix_tridiag(a,b,c);
T = diag(a) + diag(b,1) + diag(c,-1);
norm(T - hmatrix_full(P1))


P2 = hmatrix_tridiag(d,ee,f);
%P3 = hmatrix_sum(P1,P2);
%T = diag(a) + diag(b,1) + diag(c,-1)+ diag(d) + diag(ee,1) + diag(f,-1);

%norm(hmatrix_full(P3) - T)

%P3 = hmatrix_prod(P1,P2);
%T = T1 * T2;
%norm(hmatrix_full(P3) - T)

P3 = hmatrix_inv(P1);
norm(T1 * hmatrix_full(P3) - eye(n))
