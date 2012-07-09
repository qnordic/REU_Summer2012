function op = tfocs_prox( f, prox_f )
% OP = TFOCS_PROX( F, PROX_F )
%   combines F and PROX_F into the appropriate TFOCS-compatible object.
%
%   F is any function (with known proximity operator),
%   and PROX_F is its proximity operator, defined as:
%
%   PROX_F( Y, t ) = argmin_X  F(X) + 1/(2*t)*|| X - Y ||^2
%
%   To use this, please see the file PROX_L1 as an example
%
%   The basic layout of a file like PROX_L1 is as follows:
%   ( for the function F(X) = q*||X||_1 )
%
%       function op = prox_l1(q)
%       op = tfocs_prox( @f, @prox_f )
%
%         function v = f(x)
%           ... this function calculates the function f ...
%         end
%         function v = prox_f(y,t)
%           ... this function calculates the prox-function to f ...
%         end
%       end
%
%   Note: in the above template, the "end" statements are very important.
%
%
%   Also, users may wish to test their smooth function
%   with the script TEST_NONSMOOTH
%
%   See also prox_l1, test_nonsmooth


op = @fcn_impl;

function [ v, x ] = fcn_impl(x, t )
    if nargin < 1,
        error( 'Not enough arguments.' );
    end
    if nargin == 2,
        x  = prox_f(x,t);
    elseif nargout == 2,
        error( 'This function is not differentiable.' );
    end
    v = f(x);
end


end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.