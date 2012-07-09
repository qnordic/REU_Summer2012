function v = tfocs_normsq( x )

% TFOCS_NORMSQ    Squared norm. 
%    By default, TFOCS_NORMSQ(X) = TFOCS_DOT(X,X). However, certain
%    objects may have more efficient ways of computing this value.
%    If so, TFOCS_NORMSQ should be overloaded to take advantage of
%    this. However, the numerical equivalence to TFOCS_DOT(X,X) must
%    be preserved.

v = tfocs_dot( x, x );

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
