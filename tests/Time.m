% Funzione che prende in input dimensioni delle matrici, dimensioni dei blocchi finali delle h-matrici e numero di iterazioni
% e restituisce il tempo di esecuzione di test_CR e la norma del residuo
function [t,res]=Time(n,k,it)  
B=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
C=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
A=diag(rand(n,1))+diag(rand(n-1,1),-1)+diag(rand(n-1,1),1);
s=sum(A+B+C,2);
for i=1:n
	A(i,:)=A(i,:)/s(i);
	B(i,:)=B(i,:)/s(i);
	C(i,:)=C(i,:)/s(i);
end
for j=1:length(it)
	tic;
	[A0,B0,C0,G]=test_CR(A,B,C,it(j),k);
	t(j)=toc;
	G=full(G);
	res(j)=norm(B+A*G+C*G^2-G);
end
end


