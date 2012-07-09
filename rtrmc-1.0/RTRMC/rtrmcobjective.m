function [val grad hess] = rtrmcobjective(problem, U, H)
% FUNCTION [VAL GRAD] = RTRMCOBJECTIVE(PROBLEM, U)
% FUNCTION [VAL GRAD HESS] = RTRMCOBJECTIVE(PROBLEM, U, H)
% 
% Compute the value VAL and the gradient GRAD of the objective function for
% the low-rank matrix completion problem described in the structure
% PROBLEM, at the point U. If H is provided, then the Hessian HESS at U
% along the tangent vector H is computed too.
% 
% An inner mechanism prevents redundant computations, i.e., if
% RTRMCOBJECTIVE is called multiple times at the same point U, previously
% stored computations will be reused to accelerate the execution.
%
% The objective function implemented here is described in the following
% paper:
% 
% Nicolas Boumal and Pierre-Antoine Absil,
% RTRMC: A Riemannian trust-region method for low-rank matrix completion,
% Accepted at NIPS 2011, Granada (Spain), Dec. 2011.
%
% Nicolas Boumal, UCLouvain, May 19, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: rtrmc

    if nargout >= 3 && (nargin < 3 || any(size(H) ~= size(U)))
        error('The Hessian cannot be computed if H is not provided.');
    end
    
    % We will store common computations involved with recently visited U's
    % in this struct, to prevent redundant computations.
    persistent visitedUs;
    if isempty(visitedUs)
        visitedUs = struct();
    end
    
	X = problem.X;
    C = problem.C;
    m = problem.m;
    n = problem.n;
    r = problem.r;
    k = problem.k;
    I = problem.I;
    J = problem.J;
    mask = problem.mask;
	Chat = problem.Chat;
	lambda = problem.lambda;
    
    
    % Generate a field name associated to U and the problem structure
    fldname = ['z' hash([U(:); m; n; r; k; lambda; problem.id], 'md5')];
    if isfield(visitedUs, fldname)
        compumem = visitedUs.(fldname);
    else
        compumem = struct('dob', clock());
    end
    
    % Compute optimal matrix W, and the derivative of the optimal W for U
    % moving along H if the Hessian is required.
    if nargout >= 3
        [W dW compumem] = lsqfit(problem, U, H, compumem);
        compumem.W = W;
    end
    
    if ~isfield(compumem, 'W')
       [W compumem] = lsqfit(problem, U, compumem);
        compumem.W = W;
    end
    
    W = compumem.W;
    
    if ~isfield(compumem, 'UW')
        compumem.UW = spmaskmult(U, W, I, J);
    end
    UW = compumem.UW;
    
    
    % Compute objective value
    if ~isfield(compumem, 'val')
        compumem.val = .5*sum((C.*(UW-X)).^2) ...
            + .5*lambda.^2*sum(W(:).^2) ...
            - .5*lambda.^2*sum(UW.^2);
    end
    val = compumem.val;
    
    
    % Compute gradient if needed
    if nargout >= 2
		
        if ~isfield(compumem, 'sqlaWWt')
            compumem.sqlaWWt = lambda.^2*(W*W.');
        end
        sqlaWWt = compumem.sqlaWWt;
        
        
        if ~isfield(compumem, 'RU')
            compumem.RU = Chat.*(UW-X) - (lambda.^2)*X;
        end
        RU = compumem.RU;
        
        
        if ~isfield(compumem, 'grad')
            compumem.grad = multsparsefull(RU, W.', mask) + U*sqlaWWt;
        end
        grad = compumem.grad;
        
        
        % Compute Hessian along H if needed
        if nargout >= 3
		
            % since nargout >= 3, lsqfit has been called like so:
            % [W dW compumem] = lsqfit(problem, U, H, compumem);
            % and hence we know (but this is not good software practice)
            % that HW is built in compumem.
            % HW  = spmaskmult(H,  W, I, J);
            HW = compumem.HW;

            UdW = spmaskmult(U, dW, I, J);
            
            hess = multsparsefull(Chat.*(HW + UdW), W.', mask);
            hess = hess - U*(U.'*hess);
            hess = hess + multsparsefull(RU, dW.', mask);
            hess = hess + H*sqlaWWt;
            hess = hess + U*(lambda.^2*(W*dW.'));
            
        end
        
    end

    % Save the computed values in the visitedUs struct
    visitedUs.(fldname) = compumem;
    
    % If there are more than 'Ncompumem' entries in visitedUs,
    % remove the oldest ones. If Matlab runs out of memory on one of your
    % problems, you could try lowering this number. It is better if you can
    % get away with a value of at least 5 though.
    Ncompumem = 15;
    fields = fieldnames(visitedUs);
    nfields = length(fields);
    if nfields > Ncompumem
        ages = zeros(nfields, 1);
        for i = 1 : nfields
            fld = visitedUs.(fields{i});
            ages(i) = etime(clock(), fld.dob);
        end
        sortages = sort(ages, 1, 'ascend');
        maxage = sortages(Ncompumem);
        for i = 1 : nfields
            if ages(i) > maxage
                visitedUs = rmfield(visitedUs, fields{i});
            end
        end 
    end
    
end
