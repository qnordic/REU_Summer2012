function v = grProj(x, u)
% FUNCTION V = GRPROJ(X, U)
%
% Project vector U back onto the tangent space at X on the Grassmann
% manifold, i.e., V.'*X = 0.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: grExp grRetr

    v = u - x*(x.'*u);

end