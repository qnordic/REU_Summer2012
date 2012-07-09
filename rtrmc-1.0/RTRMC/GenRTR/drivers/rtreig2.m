function [x,stats] = rtreig2(A,params)
% RTREIG2  Compute the eigenvalue decomposition
%
% This computes the full eigenvalue decomposition of A by optimizing Brockett cost function on the orthogonal group 
% using the Riemannian Trust-Region with truncated CG inner solver. This method does not use the exponential 
% map for the retraction.
%
% x = rtreig(A) returns the eigenvectors x, sorted from smallest (corresponding to eigenvalue) to largest.
% [x,stats] = rtreig(A) returns x in addition to some statistics from the solver. See RTR for info.
%
% A should be a Hermitian matrix.
%
% rtreig(A,params) allows the user to specify parameters that are passed to the RTR solver.
%   params.x0        - an orthogonal matrix
%   params.Delta_bar - maximum trust-region radius
%   params.Delta0    - initial trust-region radius
%
% See also rtr, rtreig, rtresgev, rtrflat

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

   % set pointers for functions required by RTR
   n = size(A,1);
   N = diag(n:-1:1);
   fns.R    = @R;      % retraction
   fns.g    = @g;      % Riemannian metric
   fns.proj = @proj;   % projection onto tangent plane from nearby
   fns.f     = @(x)f(A,N,x);     % objective function
   fns.fgrad = @(x)grad(A,N,x);  % gradient of f
   fns.fhess = @(x,eta)H(A,N,x,eta);     % Hessian of f
   % set parameters for RTR algorithm
   params.max_inner = n*(n-1)/2;
   if ~isfield(params,'Delta_bar')
      params.Delta_bar = .5*pi*n/10;     % i don't remember where this number came from
   end
   if ~isfield(params,'Delta0'),
      params.Delta0   = params.Delta_bar/8;
   end
   if ~isfield(params,'x0')
      params.x0 = qf(randn(n,n));
   end
   [x,stats] = rtr(fns,params);

function retract = R(x,Om)
   % x is an orthogonal matrix 
   % Om is a skew-symmetric matrix representing the tangent vector  eta=x*Om
   % retract will be an orthogonal matrix
   for i=1:size(Om,1),
      Om(i,i) = Om(i,i) + 1;
   end
   retract = x*qf(Om);

function X = qf(B)
   [X,dummy] = qr(B,0);
   if rank(dummy) < size(B,2),
      error('qf error in rtrdsvd.m');
   end

function ez = g(x,Om1,Om2)
   % x is an orthogonal matrix 
   % Om1 is a skew-symmetric matrix representing the tangent vector  eta=x*Om1
   % Om2 is a skew-symmetric matrix representing the tangent vector zeta=x*Om2
   ez = 0;
   for i=1:size(Om1,1)
      ez = ez + Om1(:,i)'*Om2(:,i);
   end

function project = proj(x,Om)
   % x is an orthogonal matrix and Om is a skew-symmetric matrix representing the 
   % tangent vector eta=x*Om
   project = cbskew(Om);

function fval = f(A,N,x)
   % x is an orthogonal matrix and Om is a skew-symmetric matrix representing the 
   % tangent vector eta=x*Om
   xtax = x'*(A*x);
   fval = trace(xtax*N);

function S = cbskew(B)
   S = .5*(B - B');

function grad = grad(A,N,x)
   % x is an orthogonal matrix
   % grad will be a skew-symmetric matrix representating the tangent vector x*grad
   xtax = x'*(A*x);
   grad = brack(xtax,N);

function hess = H(A,N,x,Om)
   % x is an orthogonal matrix and Om is a skew-symmetric matrix representing the 
   % tangent vector eta=x*Om
   % hess is a skew-symmetric matrix representing the tangent vector x*hess
   xtax = x'*A*x;
   hess = 0.5*brack(brack(xtax,Om),N) + 0.5*brack(brack(N,Om), xtax);

function b = brack(A,B);
% Bracket operator
%
% Synopsis:
%  R = brack(A,B);
%
% Description:
%  Returns
%  B = A*B - B*A;

   b = A*B - B*A;
