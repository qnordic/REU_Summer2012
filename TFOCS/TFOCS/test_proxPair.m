function varargout   = test_proxPair( f, g, N, nTrials, force_gama,tol,break_on_bad )
% maxEr   = test_proxPair( f, g )
% 
%     will run tests to determine if f and g are really
%     Legendre-Fenchel conjugates.
% 
%     f and g should be in TFOCS form (see examples below)
%         which means they not only compute a function value f(x)
%         but that they can also return the proximity operator.
% 
%     The test points are assumed to be scalars by default
% 
%
% ...     = test_proxPair( f, g, N )
%     Same test, but now the test points will be vectors of length N
% 
% ...     = test_proxPair( f, g, x0 )
%     Test points will have the same size as x0
% 
% ...     = test_proxPair( f, g, [N or x0], nTrials )
%     will run "nTrials" (default: 10)
% 
% ...     = test_proxPair( f, g, [N or x0], nTrials, t )
%     will use the proximity operator:
%         prox_f(v) = argmin  f(x) + 1/(2t) *||x-v||^2
% 
%         If the function works for t = 1 but not for other t,
%         this helps you track down the bug (perhaps there is a t vs 1/t confusion
%         somewhere...)
% 
%         (default: t>0 is chosen randomly on every trial )
%
% ...   = test_proxPair( f, g, [N or x0], nTrials, t, tol )
%         sets the tolerance level for considering an inequality to be violated
%           (default is 1e-8)
%
% [x_bad, gamma_bad] = test_proxPair( ..., 'break' )
%       will return the offending point in the case that one of the the inequalities
%       has been violated. If no inequalities are violated, then the empty
%       matrices are returned.
% 
% 
% Examples of valid pairs (where q is any positive scalar, r is any scalar):
% 
% f   = prox_l1(q);
% g   = proj_linf(q);
% 
% f   = proj_l1(q);
% g   = prox_linf(q);
% 
% f   = proj_l2(q);
% g   = prox_l2(q);
% 
% f   = prox_hinge( q, r );
% g   = prox_hingeDual( q, r );
% 
% f   = proj_Rn;
% g   = proj_0;
% 
% f   = proj_Rplus;
% g   = proj_Rplus(-1); 
% 
% 
% Recall the definition of the conjugate function:
%   g(y) = f^*(y) = sup_x  <x,y> - f(x)
%   f(x)   = sup_y  <x,y> - f^*(y)
% 
% We use Moreau's Decomposition to generate identities. See Combettes and Wajs '05,
% lemma 2.10, http://www.ann.jussieu.fr/~plc/mms1.pdf
% 
% x = prox_( gamma*f )( x ) + gamma*prox_( g/gamma )( x/gamma )
% 
% and TFOCS uses (since t > 0)
% 
% prox(f)(t)(x) :=     argmin   f(v) + 1/(2t)||v-x||^2
%                =     argmin t*f(v) + 1/2  *||v-x||^2   [ N.B. This is NOT 1/t*argmin... ]
%                =     prox_( t*f )( x )
% 
% so the identity in terms of TFOCS prox is:
% 
% x = prox(f)(gamma)(x) + gamma*prox(g)(1/gamma)(x/gamma)           ( ID 1 )
%   = xf + xg
% 
% Also have the identity:
%     f(xf) + g(xg/gamma) = <xf,xg> / gamma                         ( ID 2 )
% where xf and xg defined above.
%     [ f(xf) + g(xg/gamma) >= <xf,xg> /gamma ] is the Fenchel-Young inequality.
% 
% Also have the identity:
%   ||x||^2/2 = gamma[   f^(gamma)(x) + g^(1/gamma)(x/gamma)   ]    ( ID 3 )
% where
%     f^t(x) = min f(v) + 1/(2t)||v-x||^2
% 
% These three identities are the three columns of errors output by this test function.


error(nargchk(2,7,nargin));
if nargin < 3 || isempty(N), N = 1; end
if ~isscalar(N)
    newX = @() 100*randn( size(N) );
else
    newX = @() 100*randn( N, 1 );
end

if nargin < 4 || isempty(nTrials), nTrials = 10; end
if nargin < 5 || isempty(force_gama),  force_gama = false; end
% Note: "gamma" mispelled on purpose, since Matlab's "gamma" function
%   creates problems when variables are named "gamma"
if nargin < 6 || isempty(tol), tol  = 1e-8; end
if nargin < 7, break_on_bad=[];end
if strfind( lower(break_on_bad), 'break')
    break_on_bad = true;
else
    break_on_bad = false;
end

fprintf('\n');
maxEr = 0;
FenchelYoungViolation = false;
for k = 1:nTrials
    if force_gama
        gama = double(force_gama);
    else
        gama = rand(1);
    end
    x  = newX();
    [vf,xf]     = f(x,gama);
    [vg,xg]     = g(x/gama,1/gama);
    xg  = xg*gama;
    
    er = (xf+xg) - x;
 
    % Test the scalar identity:
    %  er1 = f( xf ) + g( xg/gama ) - xf'*xg/gama;  % this can be inaccurate due to finite precision
    er2 = vf + vg - xf'*xg/gama;  % should be same as er2
 
    % And the other scalar identity:
    %    ||x||^2/2 = gama(   f^(gama)(x) + g^(1/gama)(x/gama)  )
    rhs = ( vf + 1/(2*gama)*norm(xf-x)^2 ) +  (vg + gama/2*norm(xg/gama-x/gama)^2 );
    lhs = norm(x)^2/2/gama;
    er3 = abs(rhs-lhs)/abs(lhs);
 
 
    % Another test: does the Fenchel-Young inequality hold? (this test is unlikely to find violations)
    %    for all x, y,   f(x) + g(y) >=   x'*y
    y  = newX();
    vf     = f(x);
    vg     = g(y);
    if vf + vg < x'*y
        FenchelYoungViolation = true;
    end
    
    % when both f and g are projections, xf and xg should be orthogonal
    fprintf('Random trial #%2d, errors are:\t%.2e,\t%.2e,\t%.2e\n', k, norm(er), er2, er3 );
    maxEr = max([norm(er),er2,er3,maxEr] );
    
    if break_on_bad && maxEr > tol
        varargout{1} = x;
        varargout{2} = gama;
        disp('Found a bad pair: terminating early');
        return;
    end
end


fprintf('Worst error was %.2e\n', maxEr );
if maxEr > tol
    disp('This is a BAD sign -- the functions are either not correctly implemented or there is a lot of roundoff error');
    disp('  Try running this again forcing t=1" to help find the source of error (see help file)');
else
    disp('This is a GOOD sign -- the functions are likely implemented correctly');
end
if FenchelYoungViolation
    disp('Found violation of Fenchel-Young inequality; this is BAD');
end
if break_on_bad
    varargout{1} = [];
    if nargout > 1, varargout{2} = []; end
else
    varargout{1} = maxEr;
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.