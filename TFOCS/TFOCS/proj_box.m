function op = proj_box(l,u)

%PROJ_BOX(l,u)    Projection onto the box { l <= x <= u }
%    If l or u is the empty matrix [], then the constraint is not
%    enforced (e.g. PROJ_BOX([],1) is the set { x <= 1 },
%       and PROJ_BOX(0) is the set { 0 <= x }  )
%    The parameters l and u may also be vectors or matrixes
%    of the same size as x.
%
%    OP = PROJ_BOX returns an implementation of this projection.

% warning: doesn't check to ensure l <= u

error(nargchk(1,2,nargin));
if nargin < 2, u = []; end
%op = @proj_Rplus_impl;
% bugfix, March 2011:
op = @(varargin)proj_Rplus_impl(l,u,varargin{:});

function [ v, x ] = proj_Rplus_impl( l, u, x, t )
v = 0;
switch nargin,
	case 3,
		if nargout == 2,
			error( 'This function is not differentiable.' );
        end
        if any( x < l ) || any( x > u )
            v = Inf;
        end
	case 4,
        if ~isempty(l)
            x = max( x, l );
        end
        if ~isempty(u)
            x = min( x, u );
        end
	otherwise,
		error( 'Wrong number of arguments.' );
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
