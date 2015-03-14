% Genero tre matrici tridiagonali con elementi positivi e tali che A+B+C sia stocastica
if (~exist('n', 'var'))
  n=30;
end
B=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
C=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
%C=B;
A=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
%A=diag(rand(n,1))+tril(rand(n,1)*rand(1,n))+triu(rand(n,1)*rand(1,n));
%B=diag(rand(n,1))+tril(rand(n,1)*rand(1,n))+triu(rand(n,1)*rand(1,n));
%C=diag(rand(n,1))+tril(rand(n,1)*rand(1,n))+triu(rand(n,1)*rand(1,n));
s=sum(A+B+C,2);
for i=1:n
	A(i,:)=A(i,:)/s(i);
	B(i,:)=B(i,:)/s(i);
	C(i,:)=C(i,:)/s(i);
end
