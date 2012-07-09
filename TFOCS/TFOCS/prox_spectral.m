function op = prox_spectral( q )

%PROX_SPECTRAL:    Spectral norm, i.e. max singular value.
%    OP = PROX_SPECTRAL( q ) implements the nonsmooth function
%        OP(X) = q * max(svd(X)).
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%
% This implementation uses a naive approach that does not exploit any
% a priori knowledge that X and G are low rank or sparse. Future
% implementations of TFOCS will be able to handle low-rank matrices 
% more effectively.

% Warning: this function has not been tested

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)prox_spectral_impl( q, varargin{:} );

function [ v, X ] = prox_spectral_impl( q, X, t )
if nargin < 2,
    error( 'Not enough arguments.' );
end
% if size(X,1) ~= size(X,2)
%     error('prox_spectral: variable must be a square matrix');
% end

if nargin == 3 && t > 0,
    [U,S,V] = svd( X, 'econ' );
    s = diag(S);
    tau = s(1);

    cs  = cumsum(s);
    ndx = find( cs - (1:numel(s))' .* [s(2:end);0] >= t * q, 1 );
    if ~isempty( ndx ),
        tau = ( cs(ndx) - t * q ) / ndx;
        s   = s .* ( tau ./ max( abs(s), tau ) );
    end
    X = U*diag(s)*V';
else
    if nargout == 2
        error( 'This function is not differentiable.' );
    end
    if issparse(X)
        tau = normest(X);
    else
        tau = norm(X);
    end
end
v = q * tau; 

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
