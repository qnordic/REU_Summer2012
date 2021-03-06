function op = smooth_quad( P, q, r )

%SMOOTH_QUAD   Quadratic function generation.
%   FUNC = SMOOTH_QUAD( P, q, r ) returns a function handle that implements
%
%        FUNC(X) = 0.5 * TFOCS_DOT( P * x, x ) + TFOCS_DOT( q, x ) + r.
%
%   All arguments are optional; the default values are P=I, q=0, r=0. In
%   particular, calling FUNC = SMOOTH_QUAD with no arguments yields
%   
%        FUNC(X) = 0.5 * TFOCS_NORMSQ( X ) = 0.5 * TFOCS_DOT( X, X ).
%
%   If supplied, P must be a scalar, square matrix, or symmetric linear
%   operator. Furthermore, it must be positive semidefinite (convex) or
%   negative semidefinite (concave). TFOCS does not verify operator 
%   symmetry or definiteness; that is your responsibility.
%   If P is a vector, then it assumed this is the diagonal part of P.
%   Note: when P is diagonal, this function can compute a proximity
%   operator.  If P is zero, then smooth_linear should be called instead.
%
%   If P is a scaling matrix, like the identity or a multiple of the identity
%   (say, P*x = 5*x), then specifying the scaling factor is sufficient (in
%   this example, 5). If P is empty, then P=1 is assumed.
%
%   See also smooth_linear.m

if nargin == 0,
    op = @smooth_quad_simple;
    return
end
if isa( P, 'function_handle' ),
    sz = P([],0);
    if ~isequal( sz{1}, sz{2} ),
        error( 'P must be square.' );
    end
elseif ~isnumeric( P ),
    error( 'P must be a scalar, matrix, or linear operator.' );
elseif isempty(P)
    P = 1;
elseif ndims( P ) > 2 || (~isvector(P) && size( P, 1 ) ~= size( P, 2 ) ),
    error( 'P must be a square matrix.' );
end
if nargin < 2 || isempty( q ),
    q = 0;
elseif numel(P) > 1 && numel(q) > 1 && length(q) ~= size(P,1),
    error( 'Dimension mismatch between p and q.' );
end
if nargin < 3 || isempty( r ),
    r = 0;
elseif numel(r) > 1 || ~isreal( r ),
    error( 'r must be a real scalar.' );
end
if isnumeric( P ),
    if isvector( P )
        if any(P) < 0 && any(P) < 0
            error(' P must be convex (minimization) or concave (maximization) but cannot be mixed');
        end
        P = P(:); % make it a column vector
        op = @(varargin)smooth_quad_diag_matrix( P, q, r, varargin{:} );
    else
        P = 0.5 * ( P + P' );
        op = @(varargin)smooth_quad_matrix( P, q, r, varargin{:} );
    end
else
    op = @(varargin)smooth_quad_linop( P, q, r, varargin{:} );
end

function [ v, g ] = smooth_quad_matrix( P, q, r, x, t )
switch nargin
    case 4,
        g = P * x + q;
        v = 0.5 *  tfocs_dot( x, g + q ) + r;
    case 5,
        n = length(x);
        x = (eye(n) + t*P)\ (x-t*q);
        g = x; % for this case, "g" is NOT the gradient
        v = 0.5 * tfocs_dot( P*x + 2*q, x ) + r;
end
%if nargin == 5,
    %error( 'Proximity minimization not supported by this function.' );
%end
% Note: we don't support proximity minimization, but there is
% nothing that prevents it theoretically.  You would need to know P^-1
% (and ideally calculate it quickly).  
% In particular, if P is diagonal, then it is easy.
%g = P * x + q;
%v = 0.5 *  tfocs_dot( x, g + q ) + r;



%  Jan 10, 2011: this function isn't correct for case 5
% function [ v, g ] = smooth_quad_diag_matrix( p, q, r, x, t )
% switch nargin
%     case 4,
%     case 5,
%         x = (1./(t*p+1)) .* x;
% end
% g = p .* x + q;
% v = 0.5 *  tfocs_dot( x, g + q ) + r;
function [ v, g ] = smooth_quad_diag_matrix( p, q, r, x, t )
switch nargin
    case 4,
        g = p .* x + q;
        v = 0.5 *  tfocs_dot( x, g + q ) + r;
    case 5,
        x = (1./(t*p+1)) .* (x-t*q);
        g = x; % for this case, "g" is NOT the gradient
        v = 0.5 * tfocs_dot( p.*x + 2*q, x ) + r;
end

function [ v, g ] = smooth_quad_linop( P, q, r, x, t )
if nargin == 5,
    error( 'Proximity minimization not supported by this function.' );
end
g = P( x, 1 ) + q;
v = 0.5 * tfocs_dot( x, g + q ) + r;

function [ v, x ] = smooth_quad_simple( x, t )
switch nargin,
    case 1,
    case 2,
        x = (1/(t+1)) * x;
end
v = 0.5 * tfocs_normsq( x );

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
