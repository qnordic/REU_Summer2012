% Numerical test to (sort of) validate the code for f, grad f and Hess f.
% We need to see a slope of 3 when a second order model is used.
%
% May 20, 2011  Nicolas Boumal, UCLouvain

TestRTRMC;

% pick a random point U on Grassmann
[U foo] = qr(randn(m, r), 0);

% pick a tangent vector at U
H = randn(m, r);
H = H - U*(U.'*H);
H = H / norm(H(:));

% model of f based on a truncated Taylor expansion at U along H
inner = @(A, B) trace(A.'*B);
[val grad hess] = rtrmcobjective(problem, U, H);
fmodel = [.5*inner(hess, H) inner(grad, H) val];

h = logspace(-8, 2);
fval = zeros(size(h));
for i = 1 : length(h)
    fval(i) = rtrmcobjval(problem, grExp(U, H, h(i)));
end
fmod = polyval(fmodel, h);

err = abs(fval-fmod);

loglog(h, err);
grid on;

disp(polyfit(log10(h(h >= 1e-3 & h <= 1e-2)), log10(err(h >= 1e-3 & h <= 1e-2)), 1))