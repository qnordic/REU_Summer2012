function test_dsvd(params,problem,k)
% TEST_DSVD   Test the RTRDSVD and IRTRDSVD drivers
%
% This method tests the RTRDSVD and IRTRDSVD drivers on a user-specified 
% truncated singular value problem.
% If the user does not specify the problem, then the matrix files/lenta.mtx is 
% used.
% 
% test_dsvd() uses default parameters.
%
% test_dsvd(params) allows user-specification of parameters:
%     params.x0        - initial iterate: x.U and x.V orthonormal
%     params.Delta_bar - maximum trust-region radius, default: infinity
%     params.Delta0    - initial trust-region radius, default: k*sqrt(3)
%     params.epsilon   - outer convergence tolerance, default: 1e-10
%     params.testrtr   - test rtrdsvd,  default: 1
%     params.testirtr  - test irtrdsvd, default: 1
%
% test_dsvd(params,problem) allows specification of the test matrix in a 
%                        Matrix Market file format.
%
% test_dsvd(params,problem,k) allows user-specification of the number of 
%                        basis vectors to be computed, default: 10
%
% See also rtr, irtr, irtrdsvd, rtrdsvd

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University, School of Computational Science
% Universite catholique de Louvain, Departement d'ingenierie mathematique

   if nargin < 1,
      params = [];
   end
   if nargin < 2,
      problem = 'files/lena.mtx';
   end
   fprintf('Loading %s ...\n',problem);
   A = mmread(problem);
   [m,n] = size(A);
   if nargin < 3,
      k = min(10,n);
   end

   % setup default params
   if ~isfield(params,'x0'),
      params.x0.U = orth(randn(m,k));
      params.x0.V = orth(randn(n,k));
      [U1,S1,V1] = svd(params.x0.U'*A*params.x0.V);
      params.x0.U = -params.x0.U*U1;
      params.x0.V =  params.x0.V*V1;
      clear U1 S1 v1;
   end
   if ~isfield(params,'Delta0'),
      params.Delta0 = k*sqrt(3);
   end
   if ~isfield(params,'Delta_bar'),
      params.Delta_bar = inf;
   end
   if ~isfield(params,'epsilon')
      params.epsilon = 1e-6;
   end
   if ~isfield(params,'testrtr') || params.testrtr ~= 0,
      testrtr = 1;
   else 
      testrtr = 0;
   end
   if ~isfield(params,'testirtr') || params.testirtr ~= 0,
      testirtr = 1;
   else 
      testirtr = 0;
   end

   % call the drivers
   rtr_ni = nan; irtr_ni = nan; 
   rtr_f = nan; irtr_f = nan;
   rtr_time = nan; irtr_time = nan;
   if testrtr,
       fprintf('Testing rtrdsvd...\n');
       [U ,S ,V , stats] =  rtrdsvd(A,k,params);
       rtr_ni = sum([stats.numinner]);
       rtr_f = sum(diag(S));
       rtr_time = sum([stats.time]);
   end
   if testirtr,
       fprintf('Testing irtrdsvd...\n');
       [Ui,Si,Vi,istats] = irtrdsvd(A,k,params);
       irtr_ni = sum([istats.numinner]);
       irtr_f = sum(diag(Si));
       irtr_time = sum([istats.time]);
   end

   % print some comparison
   fprintf('                             IRTR         RTR\n');
   fprintf('Num inner iterations:  %10d  %10d\n',irtr_ni,rtr_ni);
   fprintf('Final objective val :  %10.6f  %10.6f\n',irtr_f,rtr_f);
   fprintf('Total time (s)      :  %10.3f  %10.3f\n',irtr_time,rtr_time);

