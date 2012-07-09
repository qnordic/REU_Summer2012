function varargout = irtrdsvd(A,k,params)
% iRTRDSVD   Compute dominant SVD of a rectangular matrix
%
% This computes the dominant singular values and singular vectors of a rectangular matrix A
% by optimizing the function 
%   f(U,V) = trace(U'*A*V*N)
% on the manifold St(k,m) x St(k,n) using the Riemannian Trust-Region with truncated CG inner solver.
%
% The manifold St(k,m) x St(k,n) is the product of two orthogonal Stiefel manifolds:
%   St(k,n) = {X \in R^{n x k} such that X^T X = I_k}
%
% S = irtrdsvd(A,k) returns the dominant k singular values
% [U,S,V] = irtrdsvd(A,k) returns the dominant k singular vectors and values
% [U,S,V,stats] = irtrdsvd(A,k) returns in addition some statistics from the solver. See IRTR for info.
%
% irtrdsvd(A,k,params) allows the user to specify parameters that are passed to the IRTR solver.
%   params.x0        - initial iterate: x.U orthonormal and x.V orthonormal
%   params.Delta_bar - maximum trust-region radius
%   params.Delta0    - initial trust-region radius
%   params.epsilon   - Outer Convergence tolerance (absolute)
%
% See also irtr, rtreig, rtreig2, rtrflat

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University, School of Computational Science
% Universite catholique de Louvain, Departement d'ingenierie mathematique

   % set pointers for functions required by IRTR
   [m,n] = size(A);
   if nargin < 2,
      k = min(6,n);
   end
   if nargin < 3,
      params = struct([]);
   end
   if m == 0 || n == 0,
      error('First argument ''A'' to IRTRDSVD must not be empty.');
   end
   if k < 0,
      error('Second argument ''k'' to IRTRDSVD must be positive.');
   end

   N = diag([k:-1:1]);

   fns.R    = @(x,eta)R(A,x,eta);      % retraction
   fns.g    = @g;      % Riemannian metric
   fns.proj = @proj;   % projection onto tangent plane from nearby
   fns.randT = @randT;

   fns.f     = @(x)f(N,x);      % objective function
   fns.fgrad = @(x)grad(N,x);   % gradient of f
   fns.fhess = @(x,eta)H(A,N,x,eta); % Hessian of f
   fns.rho = @(x,eta)rho(A,N,x,eta);
   % set parameters for IRTR algorithm
   % use tcg solve for now
   if nargin < 3,
      params = [];
   end
   params.solver = 'tCG';
   d = (m-k)*k + k*(k-1)/2 + (n-k)*k + k*(k-1)/2;
   if ~isfield(params,'max_inner') 
      params.max_inner = d;
   else
      params.max_inner = min(params.max_inner,d);
   end
   if ~isfield(params,'Delta_bar')
      params.Delta_bar = inf;
   end
   if ~isfield(params,'Delta0'),
      params.Delta0    = k*sqrt(3);
   end
   if ~isfield(params,'x0')
      params.x0.U = qf(randn(m,k));
      params.x0.V = qf(randn(n,k));
      [U1,S1,V1] = svd(params.x0.U'*A*params.x0.V);
      params.x0.U = -params.x0.U*U1;
      params.x0.V =  params.x0.V*V1;
      clear U1 V1;
   end
   Tzero = zeros(m+n,k);
   params.x0 = R(A,params.x0,Tzero);
   clear Tzero;

   [x,stats] = irtr(fns,params);
   if nargout <= 1,
      varargout{1} = x.S;
   else
      varargout{1} = x.U;
      if nargout >= 2,
         varargout{2} = diag(x.S);
      end
      if nargout >= 3,
         % IRTR DSVD minimizes; computes (-U,V)
         % flip sign of right singular vectors
         varargout{3} = -x.V;
      end
      if nargout >= 4,
         varargout{4} = stats;
      end
   end

function X = qf(B)
   k = size(B,2);
   [X,R] = qr(B,0);
   diagR = diag(R);
   diagR(diagR==0) = 1;
   X = X*spdiags(sign(diagR),0,k,k);
   if rank(R) < size(B,2),
      error('qf error in irtrdsvd.m');
   end

function reta = R(A,x,eta)
   [m,k] = size(x.U);
   n = size(x.V,1);
   reta.U = qf(x.U + eta(1:m,1:k));
   reta.V = qf(x.V + eta(m+1:m+n,1:k));
   reta.AU = A'*reta.U;
   reta.AV = A*reta.V;
   reta.UAV = reta.U'*reta.AV;
   reta.S = svd(reta.UAV);

function ez = g(x,eta,zeta)
   ez = mytrace2(eta,zeta);

function S = mysym(B)
   S = .5*(B+B');

function Peta = proj(x,eta)
   [m,k] = size(x.U);
   n = size(x.V,1);
   Peta(1:m,1:k) = eta(1:m,1:k) - x.U*mysym(x.U'*eta(1:m,1:k));
   Peta(m+1:m+n,1:k) = eta(m+1:m+n,1:k) - x.V*mysym(x.V'*eta(m+1:m+n,1:k));

function fval = f(N,x)
   fval = trace(x.UAV*N);

function gradx = grad(N,x)
   gradx = [x.AV*N ; x.AU*N];
   gradx = proj(x,gradx);

function Heta = H(A,N,x,eta)
   [m,k] = size(x.U);
   n = size(x.V,1);
   Heta = [A*eta(m+1:m+n,1:k)*N - eta(1:m,1:k)*mysym(x.UAV*N) ; ...
           A'*eta(1:m,1:k)*N - eta(m+1:m+n,1:k)*mysym(x.UAV'*N) ];
   Heta = proj(x,Heta);

function r = randT(x)
   [m,k] = size(x.U);
   n = size(x.V,1);
   r = randn(m+n,k);
   r = r/norm(r);
   r = proj(x,r);

function r = rho(A,N,x,eta)
   fx = trace(x.UAV*N);
   xe = R(A,x,eta);
   fxe = f(N,xe);
   den = -g(x,eta,grad(N,x)) - 0.5*g(x,eta,H(A,N,x,eta));
   if abs((fx-fxe)/fx) < sqrt(eps),
       r = 1.0;
   else
       r = (fx - fxe) / den;
   end

function t = mytrace2(x,y)
   p = size(x,2);
   t = 0;
   for i=1:p,
       t = t + x(:,i)'*y(:,i);
   end

