function [x,stats] = rtrflat(params)
% RTRFLAT   Execution function for rtr.m, tailored for cost fns in R^n
%
% The use merely has to fill out the functions f, grad and H
% located at the end the the present file. 
%
% The functions R, g and proj should be disregarded. They have been
% preset for the case of a cost function in R^n.
%
% x = rtrflat(params) returns the computed optimizer.
%
% [x,stats] = rtrflat(params) returns x in addition to some
% statistics from the solver. See RTR for info. 
%
% rtrflat(params) allows the user to specify parameters that are passed to the RTR solver.
%   params.x0        - the initial iterate (required, as column vector)
%   params.Delta_bar - maximum trust-region radius
%   params.Delta0    - initial trust-region radius
%
% See also rtr

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

   % set parameters for RTR algorithm
   n = size(params.x0,1);
   params.max_inner = n*(n-1)/2;
   if ~isfield(params,'Delta_bar')
      params.Delta_bar = Inf; 
   end
   if ~isfield(params,'Delta0'),
      params.Delta0   = 1;
   end

   A = diag(1:n);   % used in my sample cost function
   mu = 1;  % likewise

   % set pointers for functions required by RTR
   fns.R    = @R;      % retraction
   fns.g    = @g;      % Riemannian metric
   fns.proj = @proj;   % projection onto tangent plane from nearby
   fns.f     = @(x)f(A,mu,x);     % objective function
   fns.fgrad = @(x)grad(A,mu,x);  % gradient of f
   fns.fhess = @(x,eta)H(A,mu,x,eta);     % Hessian of f

   % Call main RTR code:
   [x,stats] = rtr(fns,params);

function retract = R(x,eta)
   retract = x + eta;  % current iterate plus update vector

function ez = g(x,eta,zeta)
   ez = eta'*zeta;  % classical inner product in R^n

function project = proj(x,eta)
   project = eta;  % proj = identity since T_xR^n = R^n

function fval = f(A,mu,x)
   fval = x'*A*x + mu*(x'*x-1)^2;   % just a very simple example
   % (With my choice of A and mu, this function should have a
   % saddle points somewhere along each non-minor eigenvector of A,
   % and a global minimizer somewhere along the minor eigenvector
   % of A.)
   
function grad = grad(A,mu,x)
   % The gradient of f:
   grad = 2*A*x + 4*mu*(x'*x-1)*x;

function hess_eta = H(A,mu,x,eta)
   % The hessian operator of f:
   % H(x)[eta] = D gradf(x)[z].
   % 
   % This example illustrates how a Hessian, small but expensive to compute, may be stored.
   persistent hess_matrix prev_x
   if isempty(hess_matrix) || isempty(prev_x) || ~isequal(x,prev_x),
      % construct the hessian matrix
      % hess[eta] = 2*A*eta + 4*mu*(x'*x-1)*eta + 8*mu*x*(x'*eta);
      fprintf('GENERATING NEW HESSIAN\n');
      prev_x = x;
      hess_matrix = 2*A + 4*mu*(x'*x-1) + 8*mu*x*x';
   end
   % apply the hessian matrix to eta
   hess_eta = hess_matrix*eta;
