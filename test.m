A11 = rand(4); 
A22 = rand(4,4); 
a = rand(4,1); 
b = rand(4,1); 
c = rand(4,1); 
d = rand(4,1);

A = hmatrix_create(A11,A22,a,b,c,d);
x = mvm_hmatrix_avector(A, [ zeros(7,1); 1 ])
