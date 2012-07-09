function y = grExp(x, v, t)
% FUNCTION Y = GREXP(X, V, T)
%
% Exponential map on the Grassmann manifold.
%
%  Y = Exp_X(T*V), where X and Y are orthonormal matrices of the same size,
%  T is a scalar and V is a matrix of the same size as X representing a
%  tangent vector at X to the Grassmann manifold, i.e., V'*X = 0.
%  The default value for T is 1.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: grProj grRetr

    if nargin < 3 || isempty(t)
        t = 1;
    end

    r = size(x, 2);
    
    [U S V] = svds(t*v, r);

    s = diag(S);
    
    y = x*V*diag(cos(s)) + U*diag(sin(s));

end
