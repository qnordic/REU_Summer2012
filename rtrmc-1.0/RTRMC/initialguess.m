function U0 = initialguess(problem)
% FUNCTION U0 = INITIALGUESS(PROBLEM)
%
% Generates an initial guess of the column space of X, the matrix to be
% recovered, based on the observed entries of X and the mask pattern.
% This method was proposed by the OptSpace authors and is known as
% trimming. See for example:
% R.H. Keshavan and S. Oh.
% OptSpace: A gradient descent algorithm on the Grassman manifold for
% matrix completion. Arxiv preprint arXiv:0910.5260 v2, 2009
%
% Input:
%
% PROBLEM: A structure describing the low-rank matrix completion problem
%          to solve. Such a structure may be built using BUILDPROBLEM.
%
% Output:
%
% U0: An m-by-r orthonormal matrix spanning a column space that should be
%     closer to the true column space of the matrix to be recovered than a
%     simple random guess. U0 can be fed to RTRMC as initial guess to be
%     improved.
%
% Nicolas Boumal, UCLouvain, May 19, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: buildproblem rtrmc

    m = problem.m;
    n = problem.n;
    k = problem.k;
    r = problem.r;
    I = problem.I;
    J = problem.J;
    X = problem.X;
    
    % Average number of ratings per row / column
    avg_per_row = k/m;
    avg_per_col = k/n;
    
    % Number of ratings in each row / column
    n_in_rows = hist(I, 1:m);
    n_in_cols = hist(J, 1:n);

    % Select only those entries of Xm that are on rows and columns with not
    % more than twice the average number of ratings.
    % (updated after discussion with Sewoong Oh, April 2011)
    selec = (n_in_rows(I) <=  2*avg_per_row) & ...
            ...%(n_in_rows(I) >= .5*avg_per_row) & ...
            (n_in_cols(J) <=  2*avg_per_col);% & ...
            %(n_in_cols(J) >= .5*avg_per_col);
    
    % Create a sparse matrix with the appropriate sparsity structure and
    % rating values. SPARSE wants 'double' indices ... Don't know why.
    X0 = sparse(double(I(selec)), double(J(selec)), X(selec), m, n, nnz(selec));
    
    % Compute the r leading left singular vectors of X0.
    % This has a reasonable complexity.
    [U0 S0 V0] = svds(X0, r);                             %#ok<NASGU,ASGLU>
    
end
