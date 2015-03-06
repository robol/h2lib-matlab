A11 = rand(4); 
A22 = rand(4,4); 
a = rand(4,1); 
b = rand(4,1); 
c = rand(4,1); 
d = rand(4,1);

A = hmatrix_create(A11,A22,a,b,c,d);
x = mvm_hmatrix_avector(A, [ zeros(7,1); 1 ])

n = 123;
a = randn(n,1);
b = randn(n-1,1);
c = randn(n-1,1);

P = hmatrix_tridiag(a,b,c);
T = diag(a) + diag(b,1) + diag(c,-1);
TT = zeros(n);

for i = 1 : n
  TT(:,i) = mvm_hmatrix_avector(P, [ zeros(i-1,1) ; 1 ; zeros(n-i,1) ]);
end

norm(TT  - T)
