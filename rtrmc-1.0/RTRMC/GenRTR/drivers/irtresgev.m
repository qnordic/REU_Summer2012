function [V,L,stats] = irtresgev(A,B,p,params)
% IRTRESGEV   Compute extreme eigenvectors of a positive-definite Hermitian pencil
%
% This computes the space corresponding to the smallest eigenvalues of (A,B) by optimizing the Rayleigh quotient 
% on the Grassman manifold using the Implicit Riemannian Trust-Region with truncated CG inner solver.
%
% Manifold points are represented using orthonormal matrices. This is not necessary, but it simplifies
% some terms, by removing X'*B*X and inv(X'*B*X).
%
% [V,L] = irtresgev(A,B,p) returns the extreme eigenvectors of rank p.
% [V,L,stats] = irtresgev(A,B,p) returns in addition some statistics from the solver. See RTR for info.
%
% A should be a Hermitian matrix. B should be Hermitian positive-definite or empty.
%
% irtresgev(A,B,p,params) allows the user to specify parameters that are passed to the RTR solver.
%   params.x0        - initial iterate (B-orthonormal matrix)
%   params.epsilon   - Outer Convergence tolerance (absolute)
%   params.useprec   - if non-zero, irtresgev will generate a preconditioner for the problem, based on an incomplete factorization of A. This requires a positive-definite A.
%
% See also tmesgev, irtr, rtr, rtresgev

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

% This solves the extreme symmetric eigenvalue problem on the Grassmann
   % set pointers for functions required by RTR
   n = size(A,1);
   if isempty(B),
      B = speye(n);
   end
   if nargin < 4,
      params = [];
   end
   fns.R    = @(x,eta)R(A,B,x,eta);           % retraction
   fns.g    = @(x,eta,zeta)g(B,x,eta,zeta); % Riemannian metric
   fns.dist = @(x,y)dist(B,x,y);      % distance on manifold
   fns.proj = @(x,eta)proj(B,x,eta);  % projection onto tangent plane from nearby
   fns.f     = @(x)f(A,B,x);      % objective function
   fns.fgrad = @(x)grad(A,B,x);   % gradient of f
   fns.fhess = @(x,eta)H(A,B,x,eta); % Hessian of f
   fns.rho = @(x,eta)rho(A,B,x,eta);
   if p == 1,
      fns.rhosearch = @(x,eta,delta,rho_prime)rhosearch(B,x,eta,delta,rho_prime);
   end
   % set parameters for IRTR algorithm
   % use tcg solve for now
   params.solver = 'tCG';
   d = n*p - p*(p+1)/2;
   if ~isfield(params,'max_inner') 
      params.max_inner = d;
   else
      params.max_inner = min(params.max_inner,d);
   end
   if isfield(params,'userandT') && params.userandT,
      fns.randT = @(x)randT(B,x); % random vector on tangent plane
   end
   if isfield(params,'useprec') && params.useprec,
      if issparse(A),
         % perform symamd perturbation, then incomplete cholesky
         P = symamd(A);
         [R,rank] = cholinc(A(P,P),'0');
         if rank==0,
            fns.prec  = @(x,eta)precondPR(P,R,B,x,eta); % precond with inv(Pi_Bx P R' R Pi_Bx)
         else
            fprintf('Error computing incomplete Cholesky factorization of symamd(A)... IRTR/tCG will not be preconditioned.\n');
         end
      else
         % perform full cholesky
         [R,rank] = chol(A);
         if rank==0,
            fns.prec  = @(x,eta)precondR(R,B,x,eta); % precond with inv(Pi_Bx A Pi_Bx)
         else
            fprintf('Error computing Cholesky factorization of A... IRTR/tCG will not be preconditioned.\n');
         end
      end
   end
   if ~isfield(params,'x0')
      params.x0 = randn(n,p);
   end
   params.x0 = fns.R(params.x0,zeros(n,p));
   [x,stats] = irtr(fns,params);
   AA = x'*A*x;
   BB = x'*B*x;
   AA = triu(AA) + triu(AA,1)';
   BB = triu(BB) + triu(BB,1)';
   [S,L] = eig(AA,BB);
   V = x*S;

function eta = R(A,B,x,eta)
   eta = x+eta;
   xeAxe = mysym(eta'*A*eta);
   xeBxe = mysym(eta'*B*eta);
   [V,L] = eig(xeAxe,xeBxe);
   eta = eta*V;

function ez = g(B,x,eta,zeta)
   ez = mytrace2(eta,zeta);

function d = dist(B,x,y)
   x = qf(x);
   y = qf(y);
   py = y - x*(x'*y);
   d = sqrt(sum( asin(svd(py)) .^ 2 ));

function X = qf(B)
   [X,dummy] = qr(B,0);
   if rank(dummy) < size(B,2),
      error('qf error in rtrdsvd.m');
   end

function eta = proj(B,x,eta)
   Bx = B*x;
   xBBx = Bx'*Bx;
   xBBx = triu(xBBx)+triu(xBBx,1)';
   tmp = xBBx\(Bx'*eta);
   eta = eta - Bx*tmp;

function f = f(A,B,x)
   Ax = A*x;
   f = mytrace2(x,Ax);

function grad = grad(A,B,x)
   grad = proj(B,x,2*A*x);

function Heta = H(A,B,x,eta)
   Heta = proj(B,x,2*A*eta - 2*B*eta*(x'*A*x));

function eta = randT(B,x)
   eta = proj(B,x,sqrt(eps)*randn(size(x,1),size(x,2)));

function t = mytrace2(x,y)
   p = size(x,2);
   t = 0;
   for i=1:p,
       t = t + x(:,i)'*y(:,i);
   end

function S = mysym(B)
   S = .5*(B+B');

function r = rho(A,B,x,eta)
   p = size(x,2);
   if p == 1,
      r = 1/(1+eta'*B*eta);
      return;
   end
   eBe = mysym(eta'*B*eta);
   Ax = A*x;
   M = -2*eta'*Ax - eta'*(A*eta) + (eta'*B*eta)*(x'*Ax);
   r = trace( inv(eye(p)+eBe)*M ) / trace(M); 

function tau = rhosearch(B,x,eta,delta,rho_prime)
   % for p=1
   % rho(x,eta+tau*delta) 
   %   = 1 / (1 + eBe + 2*tau*eBd + tau*tau*dBd)
   % where 
   %  eBe =   eta'*B*eta
   %  eBd =   eta'*B*delta
   %  dBd = delta'*B*delta
   % 
   % must solve quadratic:
   %  0 = 1 - 1/rho_prime + eBe + 2*tau*eBd + tau*tau*dBd
   % can assume that 1/(1+eBe) > rho_prime
   %            ---> 1/rho_prime - 1 > eBe
   %            ---> 1 - 1/rho_prime + eBe < 0
   eBe = eta'*B*eta;
   eBd = eta'*B*delta;
   dBd = delta'*B*delta;
   tau = (-2*eBd + sqrt(4*eBd*eBd - 4*dBd*(1-1/rho_prime+eBe)))/(2*dBd);

function t = precondPR(P,R,B,x,s)
   % t = precondPR(P,R,B,x,s)
   % 
   % approximately solves the equation
   %   Pi_{Bx,Bx} A Pi_{Bx,Bx} t = s
   % where A is approximated by P A P \approx R' R
   %
   % with solution
   %   t = Pi_{inv(A)Bx,Bx} inv(A) s
   % It is assumed that Pi_{Bx,Bx} s = s, i.e.
   % s is perpendicular to Bx
   %
   % solve A t = R' R t = s
   % then  A u = R' R u = B x
   Bx = B*x;
   t(P,:) = R \ ( R' \  s(P,:) );
   u(P,:) = R \ ( R' \ Bx(P,:) );
   % apply Pi_{u,Bx} to t
   t = t - u * ( (Bx'*u) \ (Bx'*t) );

function t = precondR(R,B,x,s)
   % t = precondR(A,B,Bx,s)
   % solves the equation
   %   Pi_{Bx,Bx} A Pi_{Bx,Bx} t = s
   % with solution
   %   t = Pi_{inv(A)Bx,Bx} inv(A) s
   % It is assumed that Pi_{Bx,Bx} s = s, i.e.
   % s is perpendicular to Bx
   %
   % solve A t = s
   % then  A u = B x
   Bx = B*x;
   t = R\(R'\s);
   u = R\(R'\Bx);
   % apply Pi_{u,Bx} to t
   t = t - u * ( (Bx'*u) \ (Bx'*t) );

