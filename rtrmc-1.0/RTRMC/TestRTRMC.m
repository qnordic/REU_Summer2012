% Test code for the RTRMC algorithm (low-rank matrix completion)
%
% Nicolas Boumal, UCLouvain, May 19, 2011.

clear all;
close all;
clc;

cd GenRTR; ImportGenRTR; cd ..;

%% Problem instance generation

% Dimensions of the test problem
m = 1000;                                % number of rows
n = 1000;                                % number of columns
r = 10;                                  % rank
k = 5*r*(m+n-r);                         % number of known entries

% Generate an m-by-n matrix of rank r in factored form: A*B
A = randn(m, r);
B = randn(r, n);

% Pick k (or about k) entries uniformly at random
[I J k] = randmask(m, n, k);

% Compute the values of AB at these entries
X = spmaskmult(A, B, I, J);

% Define the confidence we have in each measurement X(i)
C = ones(size(X));

% Add noise if desired
noisestd = 0;
X = X + noisestd*randn(size(X));


%% Feeding the problem instance to RTRMC

% Pick a value for lambda, the regularization parameter
lambda = 1e-6;

% Build a problem structure
problem = buildproblem(I, J, X, C, m, n, r, lambda);

profile clear;
profile on;

% Compute an initial guess
initstart = tic;
U0 = initialguess(problem);
inittime = toc(initstart);

% [Optional] If we want to track the evolution of the RMSE as RTRMC
% iterates, we can do so by specifying the exact solution in factored form
% in the problem structure and asking RTRMC to compute the RMSE in the
% options structure.
problem.A = A;
problem.B = B;

% Setup the options for the RTRMC algorithm
opts.maxouter = 300;
opts.maxinner = 40;
opts.gradtol = 1e-3;
opts.verbosity = 1;
opts.order = 2;
opts.computeRMSE = true;

% Call the algorithm
[U W stats] = rtrmc(problem, opts, U0);

profile off;
profile report;


time = zeros(length(stats));
rmse = zeros(length(stats));
for i = 1 : length(stats)
    time(i) = stats(i).time;
    rmse(i) = stats(i).dist;
end
time = inittime + cumsum(time);

semilogy(time, rmse, '.-');
xlabel('Time [s]');
ylabel('RMSE');