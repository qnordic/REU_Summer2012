function op = linop_adjoint( A )

%op = LINOP_ADJOINT( A )
%    Returns a function handle to a linear operator that is the adjoint of
%    the operator supplied.

op = @(y,mode)linop_adjoint_impl( A, y, mode );

function y = linop_adjoint_impl( A, x, mode )
switch mode,
    case 0,
        y = A(x,0);
        y([1,2]) = y([2,1]);
    case 1,
        y = A(x,2);
    case 2,
        y = A(x,1);
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
