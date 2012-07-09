function [ x, out, opts ] = solver_psdComp( Xinc, opts )

% [ x, out, opts ] = solver_psdComp( Xinc, opts )
%    Solves the PSD matrix completion problem
%        minimize (1/2)*norm(X(ij)-vv).^2
%        s.t.     X p.s.d
%    where ij is a vector of indices corresponding to the known elements.
%    The nonzero values of Xinc are assumed to be the known values; that
%    is, all zero values are considered unknowns. In order to specify a
%    known zero value, replace it with a very small value; e.g., 1e-100.

% Supply default values
error(nargchk(1,2,nargin));
if nargin < 2, opts = []; end
if ~isfield( opts, 'restart' ), 
    opts.restart = 50; 
end

[n,m] = size(Xinc);
if n ~= m, error( 'Input must be square.' ); end
Xinc = tril(Xinc);
[ii,jj,vv] = find(Xinc);
ij = ii + (n-1) * jj;
linop = @(varargin)samp_op( n, ii, jj, ii+n*(jj-1), varargin{:} );
Xinc = Xinc + tril(Xinc,-1)';

% Extract the linear operators
[x,out,opts] = tfocs( smooth_quad, { linop, -vv }, proj_psd, Xinc, opts );

% TODO: update this using the new linop_subsample.m code
function y = samp_op( n, ii, jj, ij, X, mode )
switch mode
    case 0,
        y = { [n,n], [length(ii),1] };
    case 1,
        y = X(ij);
    case 2,
        y = sparse( ii, jj, X, n, n );
        y = y + tril(y,-1)';
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
