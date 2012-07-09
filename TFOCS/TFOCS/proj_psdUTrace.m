function op = proj_psdUTrace( q )

%PROJ_PSDUTRACE   Positive semidefinite cone with fixed trace.
%    OP = PROJ_PSDUTRACE( q ) returns a function that implements the
%    indicator for the cone of positive semidefinite matrices with
%    fixed trace: { X | min(eig(X+X'))>=0, trace(0.5*(X+X'))<=q
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%
% This version uses a dense eigenvalue decomposition; future versions
% of TFOCS will take advantage of low-rank and/or sparse structure.
%
% See also proj_simplex.m (the vector-analog of this function)

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
q = proj_simplex( q );
%op = @(x,t)proj_psdUTrace_q( q, x, t );
op = @(varargin)proj_psdUTrace_q( q, varargin{:} );

function [ v, X ] = proj_psdUTrace_q( eproj, X, t )

VECTORIZE   = false;
if size(X,1) ~= size(X,2)
    %error('proj_psdUTrace requires a square matrix as input');
    n = sqrt(length(X));
    X = reshape(X, n, n );
    VECTORIZE   = true;
end
v = 0;
X = 0.5*(X+X');
if nargin > 2 && t > 0,
    [V,D]=eig(X);
    [dum,D] = eproj(diag(D),t);
    tt = D > 0;
    V  = bsxfun(@times,V(:,tt),sqrt(D(tt,:))');
    X  = V * V';
    if VECTORIZE, X = X(:); end
elseif any(eig(X)<0),
    v = Inf;
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

