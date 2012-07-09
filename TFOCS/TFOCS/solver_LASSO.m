function [ x, out, opts ] = solver_LASSO( A, b, tau, x0, opts )

% [ x, out, opts ] = solver_LASSO( A, b, tau, x0, opts )
%    Solves the LASSO in the standard Tibshirani formulation,
%        minimize (1/2)*norm(A*x-b)^2
%        s.t.     norm(x,1) <= tau
%    A must be a linear operator or matrix, and b must be a vector. The
%    initial point x0 and the options structure opts are optional.
%
%   Note: many people use "LASSO" to refer to the problem:
%       minimize 1/2*norm(A*x-b)^2 + lambda*norm(x,1)
%   In TFOCS, this version is implemented in "solver_L1RLS"
%   (which stands for "L1 regularized Least Squares" )

% Supply default values
error(nargchk(3,5,nargin));
if nargin < 4, x0 = []; end
if nargin < 5, opts = []; end
if ~isfield( opts, 'restart' ), opts.restart = 100; end

% Extract the linear operators
[x,out,opts] = tfocs( smooth_quad, { A, -b }, proj_l1( tau ), x0, opts );

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

