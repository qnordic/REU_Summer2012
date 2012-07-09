function val = sqfrobnormfactors(A, B, X, Y)
% FUNCTION VAL = SQFROBNORMFACTORS(A, B, X, Y)
%
% Inputs:
%   A, X: m-by-r matrices
%   B, Y: r-by-n matrices
%
% Output:
%   val = ||AB-XY||^2_F (squared Frobenius norm of the m-by-n matrix AB-XY)
%
% Complexity:
%   O((m+n)*r^2)
%
% Accuracy:
%   val is obtained as the difference of two potentially large and very
%   close numbers, which is known to be error prone in fixed precision
%   arithmetic. Consequently, this function is not reliable if val is
%   expected to be very small compared to the squared Frobenius norms of
%   AB or XY (which are equivalent in that event).
%   When a potential accuracy drop is detected, the algorithm randomly
%   selects 50,000 entries of the matrix AB-XY, computes them and returns
%   their mean sum of squares as a proxy.
%
% Nicolas Boumal, UCLouvain, Sept. 6, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: rtrmc

    norm1 = trace((X'*X)*(Y*Y'));
    norm2 = trace((A'*A)*(B*B'));
    innerprod = trace((A'*X)*(Y*B'));
    val = norm1+norm2-2*innerprod;

    if val < 1e1 * min(eps(norm1), eps(norm2))
%         warning('RTRMC:sqfrobnormfactors', ...
%                 'The output of sqfrobnormfactors might not be reliable.');
%         disp('sqfrobnormfactors: stochastic computation :/');
        m = size(A, 1);
        n = size(B, 2);
        [I J k] = randmask(m, n, min(m*n, 50000), 1);
        AB = spmaskmult(A, B, I, J);
        XY = spmaskmult(X, Y, I, J);
        val = sum((AB-XY).^2)*m*n/k;
    end
    
end
