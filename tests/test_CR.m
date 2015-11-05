% Visualizza l'andamento del rango di quasi separabilità dell'approssimazione della matrice G, soluzione di -B+(I-A)*X-C*X^2=0
% durante l'esecuzione della cyclic reduction sulla matrice a blocchi
%
%|(I-A)   -C     0    0 0...| | G |  | B |
%|  -B  (I-A)   -C    0 0...| |G^2| =| 0 | 
%|   0    -B   (I-A) -C 0...| |G^3|  | 0 |
%|..........................| |...|  |...|
function [A0,B0,C0,G]=test_CR(A,B,C,it,bd)
n=size(A,1);
rc=Cluster(n, bd);

%a=(C-B)*ones(n,1);
%[alpha,D]=eig((A+B+C)');
%D=sum(D,2);
%[temp,i]=max(D);
%temp;
%clear temp D;
%alpha=alpha(:,i);
%alpha=alpha/sum(alpha);

B0=HMatrix('tridiagonal',rc,rc,diag(B),diag(B,-1),diag(B,1));
BB=HMatrix('tridiagonal',rc,rc,diag(B),diag(B,-1),diag(B,1));
A0=HMatrix('tridiagonal',rc,rc,ones(n,1)-diag(A),-diag(A,-1),-diag(A,1));
C0=HMatrix('tridiagonal',rc,rc,diag(C),diag(C,-1),diag(C,1));
AC=HMatrix('tridiagonal',rc,rc,ones(n,1)-diag(A),-diag(A,-1),-diag(A,1));

for i=1:it
    
    %AA = inv(A0);
	%B1=B0*(AA * B0);
	%C1=C0*(AA * C0);
	%A1=A0 - B0*(AA * C0) - C0*(AA * B0);
	%AC=AC-C0*(AA * B0);    
    % AA = inv(A0);
    % tB = AA * B0;
    % tC = AA * C0;
    t = A0 \ [ B0, C0 ];
    % t = [ tB, tC ];
    B1 = B0 * t(1);
    C1 = C0 * t(2);
    A1 = A0 - B0*t(2) - C0*t(1);
    AC = AC - C0*t(1);
    %fprintf ('A1 rank: %d\n', hmatrix_rank(A1));
	A0=A1;
	B0=B1;
	C0=C1;	
    
    hmatrix_rank(A0)

end 

	G=AC \ BB;


end
