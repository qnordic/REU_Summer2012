function op = proj_l1( q )

%PROJ_L1   The scaled 1-norm ball.
%    OP = PROJ_L1( Q ) returns an operator implementing the 
%    indicator function for the 1-norm ball of size q,
%    { X | norm( X, 1 ) <= q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, it must be a positive
%    real scalar.

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_l1_q(q, varargin{:} );

function [ v, x ] = proj_l1_q( q, x, t )
v = 0;
switch nargin,
case 2,
	if nargout == 2,
		error( 'This function is not differentiable.'  );
	elseif norm( x(:), 1 ) > q,
		v = Inf;
	end
case 3,
    s      = sort(abs(nonzeros(x)),'descend');
    cs     = cumsum(s);
    ndx    = find( cs - (1:numel(s))' .* [ s(2:end) ; 0 ] >= q, 1 );
    if ~isempty( ndx )
        thresh = ( cs(ndx) - q ) / ndx;
        x      = x .* ( 1 - thresh ./ max( abs(x), thresh ) );
    end
otherwise,
    error( 'Not enough arguments.' );
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
