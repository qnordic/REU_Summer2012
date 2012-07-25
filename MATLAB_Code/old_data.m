clear all; close all; clc;

% Loads the data from a full text file produced by a console based JAVA
% applet. The file oldDataMatrix.txt has all the 0's and 1's.

% load('oldDataMatrix.txt') ;
% A = sparse(oldDataMatrix) ;
% spy(A)

% First decomposes the matrix A into the sparse row/column form then writes
% the matrix to a text file in the sparse form.
% [i,j,val] = find(A) ;
% M = [i,j,val] ;
% fid = fopen('scirate.txt','w') ;
% fprintf(fid, '%d %d %d\n', M') ;
% fclose(fid) ;

% Loads the sprase matrix decomposition to memory and reproduces the sparse
% matrix A
M = load('scirate.txt') ;
M = spconvert(M) ;
% M = full(M) ;
% Compute the singular values of S
% [U,S,V] = svd(full(M)) ;

n1 = size(M,1) ;
n2 = size(M,2) ;
r = 10;

df = r*(n1+n2-r);
m = min(5*df,round(.99*n1*n2) );
p  = m/(n1*n2);

Omega = find(M);

data = M(Omega);

tau = 5*sqrt(n1*n2); 
delta = 1.2/p;
maxiter = 250; 
tol = 1e-4;

tic
[U,S,V,numiter,out] = SVT([n1 n2],Omega,data,tau,delta,maxiter,tol);
toc

X = U*S*V';

M = full(M) ;

% Show results
fprintf('The recovered rank is %d\n',length(diag(S)) );
fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))

Xtemp = X ;
TOL = 0.5 ;
for ii = 1:size(Xtemp,1)
    for jj = 1:size(Xtemp,2)
        if Xtemp(ii,jj) < TOL
            Xtemp(ii,jj) = 0 ;
        end
    end
end

% Xrecov = Xtemp ;
% Xrecov(find(M)) = 1 ;

% for ii = 1:size(Xtemp,1)
%     Xrecov(ii,:) = Xrecov(ii,:) - sum(Xrecov(ii,:))/length(Xrecov(ii,:)) ;
%     for jj = 1:size(Xtemp,2)
%         if Xrecov(ii,jj) < TOL
%             Xrecov(ii,jj) = 0 ;
%         end
%     end
% end

xlswrite('SVT_sol.xlsx', Xrecov)