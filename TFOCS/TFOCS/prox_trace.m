function op = prox_trace( q, LARGESCALE )

%PROX_TRACE     Nuclear norm, for positive semidefinite matrices. Equivalent to trace.
%    OP = PROX_TRACE( q ) implements the nonsmooth function
%        OP(X) = q * sum(svd(X)) = q*tr(X) ( X >= 0 assumed )
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%
%    OP = PROX_TRACE( q, LARGESCALE )
%       uses a Lanczos-based Eigenvalue decomposition if LARGESCALE == true,
%       otherwise it uses a dense matrix Eigenvalue decomposition
%
%    CALLS = PROX_TRACE( 'reset' )
%       resets the internal counter and returns the number of function
%       calls
%
% This implementation uses a naive approach that does not exploit any
% a priori knowledge that X and G may be low rank (plus sparse). Future
% implementations of TFOCS will be able to handle low-rank matrices 
% more effectively.

if nargin == 1 && strcmpi(q,'reset')
    op = prox_trace_impl;
    return;
end

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
if nargin < 2, LARGESCALE = []; end

% clear the persistent values:
prox_trace_impl();

op = @(varargin)prox_trace_impl( q, varargin{:} );


function [ v, X ] = prox_trace_impl( q, X, t )
persistent oldRank
persistent nCalls
persistent V
if nargin == 0, oldRank = []; v = nCalls; nCalls = []; V=[]; return; end
if isempty(nCalls), nCalls = 0; end

if nargin > 2 && t > 0,
    
    if ~isempty(LARGESCALE) % inherited from parent
        largescale = LARGESCALE;
    else
        largescale = ( numel(X) > 100^2 ) && issparse(X);
    end
    tau = q*t;
    nCalls = nCalls + 1;
    
    if ~largescale
        [V,D]   = eig(full((X+X')/2));
    else
        %% Experimental code:
        %if any(any(isinf(X))), disp('found Inf!'); keyboard; end
        %
        %opts = [];
        %opts.lambda = tau/1.1; % be conservative!
        %opts.iters  = 4;
        %[V,D] = eigs_approximate( X, V, opts );
        %fprintf('Size of V is %d\n', size(V,2) );
%%         er=norm(V*D*V' - X,'fro')/norm(X,'fro');
%%         fprintf('Size of V is %d and approx error is %.2e\n', size(V,2), er);
%%         if er > 1.99, disp('BREAK!'); keyboard; end
        %
    %end
    %if false
        
        % Guess which eigenvalue value will have value near tau:
        [M,N] = size(X);
        X = (X+X')/2;
        if isempty(oldRank), K = 10;
        else, K = oldRank + 2;
        end
        
        ok = false;
        opts = [];
        opts.tol = 1e-10; 
        if isreal(X)
            opts.issym = true;
        end
        while ~ok
            K = min( [K,M,N] );
            [V,D] = eigs( X, K, 'LM', opts );
            ok = (min(diag(D)) < tau) || ( K == min(M,N) );
            if ok, break; end
            K = 2*K;
            if K > 10
                opts.tol = 1e-6;
            end
            if K > 40
                opts.tol = 1e-4;
            end
            if K > 100
                opts.tol = 1e-3;
            end
            if K > min(M,N)/2
                [V,D]   = eig(full((X+X')/2));
                ok = true;
            end
        end
        oldRank = length(find(diag(D) > tau));
    end
    s  = diag(D) - tau;
    tt = s > 0;
    s  = s(tt,:);
    
    if isempty(s),
        X = tfocs_zeros(X);
    else
        X = V(:,tt) * bsxfun( @times, s, V(:,tt)' );
    end
    v = q * sum(s);

else
    v = q* trace( X + X' )/2;
end

end

end % end of main function

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
