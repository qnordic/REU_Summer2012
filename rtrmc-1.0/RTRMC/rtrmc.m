function [U W stats] = rtrmc(problem, opts, U0)
% FUNCTION [U W STATS] = RTRMC(PROBLEM, OPTS, U0)
%
% A Riemannian Trust-Region (RTR) method for Matrix Completion.
%
% Inputs:
%
% PROBLEM: A structure describing the low-rank matrix completion problem to
%          solve. Such a structure can be constructed using BUILDPROBLEM.
%
% OPTS:    A structure holding a few parameters about the optimization
%          algorithm to be used. The parameters are:
%           *maxouter: maximum number of outer iteration of the RTR method
%           *maxinner: maximum number of inner iterations of the RTR method
%           *gradtol: tolerance on the gradient norm (stopping criterion)
%           *verbosity: 0 for silent mode, 1 to display iteration information
%           *order: 1 for RTRMC 1 (w/o Hessian) and 2 for RTRMC 2 (w/ Hessian)
%           *computeRMSE: boolean. If the problem structure contains matrices
%            A and B such that the true matrix to recover is A*B and if
%            this flag is set to true, at each iteration the Root Mean
%            Squared Error between A*B and U*W is computed and put in the
%            STATS structure in the field DIST.
%
%            All parameters are optional.
%
% U0:      Orthonormal matrix of size m-by-r, initial guess for the column
%          space spanned by the matrix to recover. An initial guess may be
%          computed using INITIALGUESS.
% 
% Outputs:
% 
% U, W: U is an m-by-r orthonormal matrix, W is an r-by-n matrix.
%       The reconstructed matrix is U*W.
% 
% STATS: Structure containing information about the optimization iterations.
% 
% The algorithm implemented here is described in the following paper:
% 
% Nicolas Boumal and Pierre-Antoine Absil,
% RTRMC: A Riemannian trust-region method for low-rank matrix completion,
% Accepted at NIPS 2011, Granada (Spain), Dec. 2011.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: buildproblem initialguess


    if nargin < 3 || isempty(U0)
        U0 = initialguess(problem);
    end
    

    % Default option values
    if ~isfield(opts, 'maxouter'),      opts.maxouter = 300; end
    if ~isfield(opts, 'maxinner'),      opts.maxinner = 40; end
    if ~isfield(opts, 'gradtol'),       opts.gradtol = 1e-5; end
    if ~isfield(opts, 'verbosity'),     opts.verbosity = 1; end
    if ~isfield(opts, 'order'),         opts.order = 2; end
    if ~isfield(opts, 'computeRMSE'),   opts.computeRMSE = false; end
    
    
    % Check some of the option values
    if opts.order ~= 1 && opts.order ~= 2
        warning('RTRMC:order', ...
          ['The RTRMC ''order'' option must be set to 1 (no Hessian) ' ...
           'or 2 (with Hessian).']);
        opts.order = 2;
    end
    

    % retraction
    fns.R     = @grRetr;

    % Riemannian metric
    fns.g     = @(U, H1, H2) trace(H1.'*H2);

    % Projection onto tangent plane from nearby
    fns.proj  = @(U, H) H - U*(U.'*H);

    % If the problem structure contains a pair (A, B) such that the true
    % matrix we are looking for is AB, RTRMC can compute the total RMSE
    % at each iteration efficiently. 
    if opts.computeRMSE && isfield(problem, 'A') && isfield(problem, 'B')
        fns.dist = @computeRMSE;
        params.xsol = problem;
    end

    % Returns a very small random tangent vector in T_x M
    fns.randT = @(U) fns.proj(U, 1e-10*randn(size(U)));

    % Objective function
    fns.f     = @(U) rtrmcobjval(problem, U);

    % Gradient of f
    fns.fgrad = @(U) rtrmcobjgrad(problem, U);

    % Hessian of f
    if opts.order == 2
        fns.fhess = @(U, H) rtrmcobjhess(problem, U, H);
    elseif opts.order == 1
        % If the Hessian is not available, 'approximate' the Hessian
        % by the identity map.
        fns.fhess = @(U, H) H;
    end

    params.x0 = U0;
    params.Delta_bar = .5;
    params.Delta0   = params.Delta_bar/8;
    params.max_inner = opts.maxinner;
    params.max_outer = opts.maxouter;
    params.useRand = 1;
    params.testgh = 0;
    params.debug = 0;
    params.verbosity = opts.verbosity;
    params.epsilon = opts.gradtol;

    [U stats] = rtr(fns, params);
    
    W = lsqfit(problem, U);

end

function rmse = computeRMSE(U, problem)

    if isfield(problem, 'A') && isfield(problem, 'B')
        
        W = lsqfit(problem, U);
        A = problem.A;
        B = problem.B;
        m = problem.m;
        n = problem.n;
        rmse = sqrt(sqfrobnormfactors(U, W, A, B)/(m*n));
        
    else
        rmse = NaN;
    end
    
end
