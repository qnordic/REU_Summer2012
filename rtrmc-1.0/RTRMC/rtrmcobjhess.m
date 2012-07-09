function hess = rtrmcobjhess(problem, U, H)
% Helper function for RTRMC
    
    [val grad hess] = rtrmcobjective(problem, U, H); %#ok<ASGLU>

end
