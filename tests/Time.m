% Funzione che prende in input dimensioni delle matrici, dimensioni dei blocchi finali delle h-matrici e numero di iterazioni
% e restituisce il tempo di esecuzione di test_CR e la norma del residuo
function [t,res]=Time(n,p,q,k,it)  
[A,B,C]=bandatirand(n,p,q);
for j=1:length(it)
	tic;
	[A0,B0,C0,G]=test_CR(A,B,C,it(j),k);
	t(j)=toc;
	G=full(G);
	res(j)=norm(B+A*G+C*G^2-G);
end
end


