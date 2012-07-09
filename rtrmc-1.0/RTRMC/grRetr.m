function y = grRetr(x, v, t)
% FUNCTION Y = GRRETR(X, V, T)
%
% Retraction map on the Grassmann manifold.
%
%  Y = Retr_X(T*V), where X and Y are orthonormal matrices of the same
%  size, T is a scalar and V is a matrix of the same size as X representing
%  a tangent vector at X to the Grassmann manifold, i.e., V'*X = 0.
%  The default value for T is 1.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: grProj grExp

    if nargin < 3 || isempty(t)
        t = 1;
    end
    
    [y R] = qr(x+t*v, 0);         %#ok<NASGU>

end
