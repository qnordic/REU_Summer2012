function test_esgev(params,problem,p)
% TEST_ESGEV   Test the RTRESGEV, IRTRESGEV and TMESGEV drivers
%
% This method tests the RTRESGEV, IRTRESGEV and TMESGEV drivers on a user-specified 
% generalized eigenvalue problem.
% If the user does not specify the problem, then the Matrix Market problem 
% BCSST_08 is used.
% 
% test_esgev() uses default parameters.
%
% test_esgev(params) allows user-specification of parameters:
%     params.x0        - initial iterate: x.U and x.V orthonormal
%     params.rho_prime - trust-region acceptance parameters, default: .1
%     params.useprec   - use a preconditioner (see rtresgev), default: 0
%     params.testtm    - test tmesgev,   default: 1
%     params.testrtr   - test rtresgev,  default: 1
%     params.testirtr  - test irtresgev, default: 1
%
% test_esgev(params,problem) allows specification of the test matrix in a 
%                        Matrix Market file format.
%                        problem or problem{1} should specify filename for A
%                        problem{2} should specify filename for B
%                        default: problem = {'files/bcsstk08.mtx',
%                                            'files/bcsstm08.mtx'}
%
% test_esgev(params,problem,p) allows user-specification of the number of 
%                        basis vectors to be computed, default: 6
%
% See also rtr, irtr, irtresgev, rtresgev, tmesgev

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University, School of Computational Science
% Universite catholique de Louvain, Departement d'ingenierie mathematique

   if nargin < 1,
      params = [];
   end
   if nargin < 2,
      problem{1} = 'files/bcsstk08.mtx';
      problem{2} = 'files/bcsstm08.mtx';
   end
   if iscell(problem),
      probA = problem{1};
      if length(probA) > 1,
         probB = problem{2};
      else 
         probB = [];
      end
   end
   fprintf('Loading ''%s'' for A...\n',probA);
   K = mmread(probA);
   M = [];
   if ~isempty(probB),
      fprintf('Loading ''%s'' for B...\n',probB);
      M = mmread(probB);
   end
   [m,n] = size(K);
   if nargin < 3,
      p = min(6,n);
   end

   % initialize both with same starting subspace
   if ~isfield(params,'x0'),
      params.x0 = randn(n,p);
   end
   if ~isfield(params,'rho_prime'),
      params.rho_prime = .1;
   end
   if ~isfield(params,'testtm') || params.testtm ~= 0,
      testtm = 1;
   else 
      testtm = 0;
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
   rtr_ni = nan; irtr_ni = nan; tm_ni = nan;
   rtr_f = nan; irtr_f = nan; tm_f = nan;
   rtr_time = nan; irtr_time = nan; tm_time = nan;
   if testrtr,
       fprintf('Testing RTR-ESGEV...\n');
       [V ,S , stats] =  rtresgev(K,M,p,params);
       rtr_ni = sum([stats.numinner]);
       rtr_f = sum(diag(S));
       rtr_time = sum([stats.time]);
   end
   if testirtr,
       fprintf('Testing IRTR-ESGEV...\n');
       [Vi,Si,istats] = irtresgev(K,M,p,params);
       irtr_ni = sum([istats.numinner]);
       irtr_f = sum(diag(Si));
       irtr_time = sum([istats.time]);
   end
   if testtm,
       fprintf('Testing TM-ESGEV...\n');
       [Vt,St,tmstats] =   tmesgev(K,M,p,params);
       tm_ni = sum([tmstats.numinner]);
       tm_f = sum(diag(St));
       tm_time = sum([tmstats.time]);
   end
   fprintf('                             IRTR         RTR          TM\n');
   fprintf('Num inner iterations:  %10d  %10d  %10d\n',irtr_ni,rtr_ni,tm_ni);
   fprintf('Final objective val :  %10.6f  %10.6f  %10.6f\n',irtr_f,rtr_f,tm_f);
   fprintf('Total time (s)      :  %10.3f  %10.3f  %10.3f\n',irtr_time,rtr_time,tm_time);
