function op = proj_simplex( q )

%PROJ_SIMPLEX	Projection onto the simplex.
%    OP = PROJ_SIMPLEX( Q ) returns an nonsmooth function that
%    represents the scaled simplex { x | x >= 0, sum(x) <= q }.
%    Q is optional; if not supplied, it defaults to 1. If it  is
%    supplied, it must be a real positive scalar.
%
%   See also proj_psdUTrace.m (the matrix-analog of this function)

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_simplex_q( q,varargin{:} );

function [ v, x ] = proj_simplex_q( q, x, t )
v = 0;
if nargin > 2 && t > 0,
	if any( x < 0 ) || sum( x ) > q,
        s     = sort( x, 'descend' );
        cs    = ( cumsum(s) - q ) ./ ( 1 : numel(s) )';
        ndx   = nnz( s > cs );
        x     = max( x - cs(ndx), 0 );
	end
elseif any( x < 0 ) || abs( sum(x) / q - 1 ) > sqrt(numel(x)) * eps,
    v = Inf;
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
