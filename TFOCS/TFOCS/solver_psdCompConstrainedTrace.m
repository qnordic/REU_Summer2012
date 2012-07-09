function [ x, out, opts ] = solver_psdCompConstrainedTrace( Xinc, tr, opts )

% [ x, out, opts ] = solver_psdCompConstrainedTrace( Xinc, tr, opts )
%    Solves the PSD matrix completion problem
%        minimize (1/2)*norm(X(ij)-vv).^2
%        s.t.     X p.s.d and trace(X) = tr
%    where ij is a vector of indices corresponding to the known elements.
%    The nonzero values of Xinc are assumed to be the known values; that
%    is, all zero values are considered unknowns. In order to specify a
%    known zero value, replace it with a very small value; e.g., 1e-100.
%
%   See also solver_psdComp

% Supply default values
error(nargchk(1,3,nargin));
if nargin < 3, opts = []; end
if ~isfield( opts, 'restart' ), 
    opts.restart = 50; 
end
if nargin < 2 || isempty(tr)
    tr = 1;
end

[n,m] = size(Xinc);
if n ~= m, error( 'Input must be square.' ); end
Xinc = tril(Xinc);
[ii,jj,vv] = find(Xinc);
ij = ii + (n-1) * jj;
linop = @(varargin)samp_op( n, ii, jj, ii+n*(jj-1), varargin{:} );
Xinc = Xinc + tril(Xinc,-1)';

% Extract the linear operators
[x,out,opts] = tfocs( smooth_quad, { linop, -vv }, proj_psdUTrace(tr), Xinc, opts );

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
