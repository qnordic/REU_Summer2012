function op = smooth_entropy()
op = @smooth_entropy_impl;

function [ v, g ] = smooth_entropy_impl( x )
if any( x < 0 ),
    v = -Inf;
    if nargout > 1,
        g = NaN * ones(size(x));
    end
else
    logx = log(max(x,realmin));
    v = - tfocs_dot( x, logx );
    if nargout > 1,
        g = - logx - 1;
    end
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
