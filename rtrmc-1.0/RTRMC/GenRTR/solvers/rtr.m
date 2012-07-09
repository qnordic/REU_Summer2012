function varargout = rtr(fns,params)
% RTR   Riemannian Trust-Region (with tCG inner solve)
%
% This performs a retraction-based, trust-region minimization of an objective function on a 
% Riemannian manifold.
%
% x = rtr(fns,params) returns the minimizer (a point on the manifold)
% [x,stats] = rtr(fns,params) also returns some statistics from the algorithm
%
% The struct fns contains the function handles necessary to perform the optimization:
%  fns.R(x,eta) : retract tangent vector eta at T_x M to the manifold
%  fns.g(x,eta,zeta) : Riemannian metric at T_x M 
%  fns.proj(x,eta) : Project eta from nearby T_x M to T_x M (used in iteration to combat round-off)
%                    Also used in debugging to check tangentiality
%                    Can be set to the identity if this is not a concern.
%  fns.f(x) : Compute the objective function at x
%  fns.fgrad(x) : Compute the gradient of f at x, returning a tangent vector
%  fns.fhess(x,eta) : Apply the Hessian of f, hess_x: T_x M -> T_x M
% Optionally, three other methods may be specified:
%  fns.prec(x,eta) : apply an s.p.d. preconditioner for the inner iteration, to eta, a
%                    tangent vector in T_x M
%  fns.dist(x,y) : returns the distance on the manifold between points x and y
%  fns.randT(x) : returns a very small random tangent vector in T_x M
% More info:
%  fns.dist is used to measure the distance from solution and is only referenced if params.xsol was specified.
%  fns.randT is used to initialize the trust-region subproblem with a random vector, to encourage escape from an
%      exact critical point.
%  fns.prec is a preconditioner for the model minimization, a positive-definite (undr fns.g)
%      mapping from T_x M to T_x M
%
% Parameters to the solver include:
% Required parameters:
%   params.x0        - An initial iterate. Because RTR does not know what the manifold is, this is mandatory.
%   params.Delta_bar - Maximum trust-region radius. Because RTR does not know anything about the manifold, this is mandatory.
%   params.Delta0    - Initial trust-region radius. Because RTR does not know anything about the manifold, this is mandatory.
%   params.xsol      - For testing distance from solution using fns.dist
%   params.verbosity - Level of verbosity: 0 is silent, 1 has one-line info, 2 has detailed info
%   params.debug     - Debugging tests (0 is silent, 1 has some, 2 has a lot, 3 is more than you want to know)
%   params.min_outer - Minimum number of outer iterations (default: 0); used only with randomization
%   params.max_outer - Maximum number of outer iterations (default: 100)
%   params.min_inner - Minimum number of inner iterations (default: 0). Only in effect if using randomized initial tangent vectors.
%   params.max_inner - Maximum number of inner iterations (default: inf). Recommended: dimension of manifold.
%   params.epsilon   - Outer Convergence tolerance (absolute)
%   params.kappa     - Inner kappa convergence tolerance
%   params.theta     - Inner theta convergence tolerance
%   params.rho_prime - Accept/reject ratio
%   params.testgh    - perform some simple numerical testing of gradient and Hessian
%   params.useRand   - initial the trust-region solve with a random tangent vector
%
% The stats output is a struct array, with fields:
%  k  - the outer iteration number for the stats
%  ng - the norm of the gradient, sqrt(g(x,gradfx,gradfx))
%  fx - the current value under the objective function
%  rho - the performance ratio for the iterate
%  time - the wallclock time for the outer iteration
%  accepted - whether the proposed iterate was accepted or not
%  numinner - the number of inner iterations used to compute the next iterate
%  Delta - the trust-region radius at the outer iteration
%  dist     - the distance from the solution
%  cauchy   - whether the cauchy point was used in place of an updated computed from a random initial tangent vector
%
% See also irtr

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

% Modification history
% Version 0.221 - CGB
%   Using preconditioner-based trust-region radius
%   Disable randomization when using preconditioner
% Version 0.2 - CGB
%   Added more commenting
%   Added preconditioning for inner iteration (fns.prec)
% Version 0.1 - CGB - Sat Nov 04 2006
%   Started from RTR03(Baker), from tCG07(Absil)

tcg_stop_reason = {'negative curvature',...
                   'exceeded trust region',...
                   'reached target residual-kappa',...
                   'reached target residual-theta',...
                   'dimension exceeded'};
lan_stop_reason = {'reached residual (interior)',...
                   'reached residual (exterior)',...
                   'breakdown - reduced T',...
                   'dimension exceeded'};

if nargin < 2,
   error('Invalid arguments to rtr. See ''help rtr''');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check functions handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if   ~(isfield(fns,'R')     && isa(fcnchk(fns.R)    ,'function_handle')) ...
  || ~(isfield(fns,'g')     && isa(fcnchk(fns.g)    ,'function_handle')) ...
  || ~(isfield(fns,'proj')  && isa(fcnchk(fns.proj) ,'function_handle')) ...
  || ~(isfield(fns,'f')     && isa(fcnchk(fns.f)    ,'function_handle')) ...
  || ~(isfield(fns,'fgrad') && isa(fcnchk(fns.fgrad),'function_handle')) ...
  || ~(isfield(fns,'fhess') && isa(fcnchk(fns.fhess),'function_handle')),
  error('Invalid arguments to rtr. See ''help rtr''');
end
fns.havedist = 0;
fns.haverand = 0;
fns.haveprec = 0;
if isfield(fns,'dist') && isa(fcnchk(fns.dist),'function_handle')
    fns.havedist = 1;
end
if isfield(fns,'randT') && isa(fcnchk(fns.randT),'function_handle')
    fns.haverand = 1;
else
    fns.haverand = 0;
end
if isfield(fns,'prec') && isa(fcnchk(fns.prec),'function_handle'),
    fns.haveprec = 1;
else
    fns.haveprec = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbosity =   1;
debug     =   0;
min_inner =   0;
max_inner = inf;
min_outer =   3;
max_outer = 100;
epsilon   =   1e-6;
kappa     =   0.1;
theta     =   1.0;
rho_prime =   0.1;
kappa_easy = .001;
testgh    = 0;
useRand  = 0;
%Delta_bar =  user must specify
%Delta0    =  user must specify
%x0        =  user must specify

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[verbosity,params] = get_int(params,'verbosity','Level of verbosity'                ,verbosity);
[debug    ,params] = get_int(params,'debug'    ,'Debugging tests'                   ,debug);
[min_outer,params] = get_int(params,'min_outer','Minimum number of outer iterations',min_outer,0);
[max_outer,params] = get_int(params,'max_outer','Maximum number of outer iterations',max_outer,0);
[min_inner,params] = get_int(params,'min_inner','Minimum number of inner iterations',min_inner,0);
[max_inner,params] = get_int(params,'max_inner','Maximum number of inner iterations',max_inner,0);
[testgh   ,params] = get_int(params,'testgh'   ,'Test gradient and Hessian'         ,testgh);
[useRand  ,params]  = get_int(params,'useRand'  ,'Random TR subproblem init'        ,useRand);
[x0       ,params] = get_iterate(params,'x0'   ,'Initial iterate');
[xsol     ,params] = get_iterate(params,'xsol' ,'Solution for testing',[]);
if isempty(xsol),
    % must have fns.dist and xsol to make fns.dist() comparisons
    fns.havedist = 0;
end
[epsilon  ,params] = get_float(params,'epsilon'  ,'Outer Convergence tolerance'      ,epsilon,0);
[rho_prime,params] = get_float(params,'rho_prime','Accept/reject parameter'          ,rho_prime,0,1);
[Delta_bar,params] = get_float(params,'Delta_bar','Maximum trust-region radius');
[Delta0   ,params] = get_float(params,'Delta0'   ,'Initial trust-region radius');
[kappa    ,params] = get_float(params,'kappa'    ,'Inner kappa convergence tolerance',kappa,0,1);
[theta    ,params] = get_float(params,'theta'    ,'Inner theta convergence tolerance',theta,0);
[kappa_easy,params] = get_float(params,'kappa_easy','Inner convergence tolerance'    ,kappa_easy,0,1);

if testgh && fns.haverand,
   TestGH(x0,fns);
   varargout{1} = x0;
   if nargout > 1,
      varargout{2} = struct([]);
   end
   return;
end

if fns.haveprec && useRand && fns.haverand,
   fprintf('rtr/tCG: Cannot perform preconditioned tCG with random initial vector.\n');
   fprintf('rtr/tCG: Disabling randomization.\n');
   useRand = 0;
end

if ~fns.haverand,
   useRand = 0;
end

% ***** Initializations *****

% initialize counters/sentinals
% allocate storage for dist, counters
k          = 0;  % counter for outer (TR) iteration.
stop_outer = 0;  % stopping criterion for TR.
total_time = 0;     % total time spent in outer loop


% initialize solution and companion measures: 
% fgrad(x)
% f(x)
tic
x = x0;
fx = fns.f(x);
fgradx = fns.fgrad(x);
norm_grad = sqrt(fns.g(x,fgradx,fgradx));
this_time = toc;

% initialize trust-region radius
Delta = Delta0;

% ** Record data:
curstat.k        = k;
curstat.ng       = norm_grad;
curstat.fx       = fx;
curstat.rho      = inf;
curstat.rhonum      = 0;
curstat.rhoden      = 0;
curstat.time = this_time;
curstat.accepted = true;
curstat.numinner = 0;
curstat.Delta    = Delta;
if fns.havedist,
    curstat.dist = fns.dist(x,xsol);
end
if useRand,
    curstat.cauchy = false;
end
stats(1) = curstat;

% ** Display:
if (verbosity == 1)
   fprintf(['%3s %3s      %5s                %5s     ',...
            'f: %e   |grad|: %e\n'],...
           '   ','   ','     ','     ',fx,norm_grad);
elseif (verbosity > 1)
   fprintf('************************************************************************\n');
   fprintf(['%3s %3s    k: %5s     num_inner: %5s     %s\n'],...
           '','','______','______','');
   fprintf('       f(x) : %e       |grad| : %e\n',fx,norm_grad);
   fprintf('      Delta : %f\n',Delta);
   if fns.havedist,
      fprintf('       dist : %e\n',stats(end).dist);
   end
   fprintf('       Time : %f\n',this_time);
end


% **********************
% ** Start of TR loop **
% **********************
while stop_outer==0,

   % start clock for this outer iteration
   tic

   % update counter
   k = k+1;

   if (verbosity > 1) || (debug > 0),
      fprintf('************************************************************************\n');
   end

   % *************************
   % ** Begin TR Subproblem **
   % *************************
  
   % determine eta0
   if useRand,
      % random vector in T_x M
      eta = fns.randT(x);
      % must be inside trust-region
      while fns.g(x,eta,eta) > Delta^2,
         % eta = eta * sqrt(sqrt(eps));
         eta = tgtveclincomb(sqrt(sqrt(eps)), eta);
      end
   else
      % without randT, 0*fgradx is the only way that we 
      % know how to create a tangent vector
      % eta = 0*fgradx;
      eta = tgtveclincomb(0, fgradx);
   end

   % solve TR subproblem
   [eta,numit,stop_inner] = tCG(fns,x,fgradx,eta,Delta,theta,kappa,min_inner,max_inner,useRand,debug);
   srstr = tcg_stop_reason{stop_inner};
   norm_eta = sqrt(fns.g(x,eta,eta));
   if debug > 0,
      testangle = fns.g(x,eta,fgradx) / (norm_eta*norm_grad);
   end

   % if using randomized approach, compare result with the
   % cauchy point
   % convergence proofs assume that we achieve at least the
   %   reduction of the cauchy point
   if useRand,
      used_cauchy = false;
      % check the curvature
      Hg = fns.fhess(x,fgradx);
      g_Hg = fns.g(x,fgradx,Hg);
      if g_Hg <= 0, 
         tau_c = 1;
      else
         tau_c = min( norm_grad^3/(Delta*g_Hg) , 1);
      end
      % and gen the cauchy point
      %eta_c = -tau_c * Delta / norm_grad * fgradx;
      eta_c = tgtveclincomb(-tau_c*Delta/norm_grad, fgradx);
      
      % now that we have computed the cauchy point in addition
      % to the returned eta, we might as well keep the better 
      % of them
      Heta = fns.fhess(x,eta);
      Heta_c = fns.fhess(x,eta_c);
      mdle  = fx + fns.g(x,fgradx,eta  ) + 1/2*fns.g(x,Heta  ,eta  );
      mdlec = fx + fns.g(x,fgradx,eta_c) + 1/2*fns.g(x,Heta_c,eta_c);
      if ( mdle > mdlec ) 
         eta = eta_c;
         norm_eta = sqrt(fns.g(x,eta,eta));
         used_cauchy = true;
      end 
   else
     Heta = fns.fhess(x,eta);
   end 

   % compute the retraction of the proposal
   x_prop  = fns.R(x,eta);

   % compute function value of the proposal
   fx_prop = fns.f(x_prop);

   % do we accept the proposed solution or not?
   % compute the Hessian at the proposal
   Heta = fns.fhess(x,eta);

   % check the performance of the quadratic model
   rhonum = fx-fx_prop;
   rhoden = -fns.g(x,fgradx,eta) - 0.5*fns.g(x,Heta,eta);
   if debug > 0,
      if rhoden < 0,
         fprintf('Error! no model decrease!\n');
         keyboard;
      end
   end
   if debug > 0,
      fprintf('DBG:     rhonum : %e\n',rhonum);
      fprintf('DBG:     rhoden : %e\n',rhoden);
   end
   rho =   rhonum  / rhoden;
   if debug > 1,
      m = @(x,eta) fns.f(x) + fns.g(x,fns.fgrad(x),eta) + .5*fns.g(x,fns.fhess(x,eta),eta);
      actrho = (fns.f(x) - fns.f(x_prop)) / (m(x,tgtveclincomb(0, eta)) - m(x,eta));
      fprintf('DBG:   new f(x) : %e\n',fx_prop);
      fprintf('DBG: actual rho : %e\n',actrho);
   end
   % HEURISTIC WARNING:
   % if abs(model change) is relatively zero, we are probably near a critical
   % point. set rho to 1.
   if abs(rhonum/fx) < sqrt(eps),
      small_rhonum = rhonum;
      rho = 1;
   else 
      small_rhonum = 0;
   end

   % choose new TR radius based on performance
   trstr = '   ';
   if rho < 1/4
      trstr = 'TR-';
      %Delta = 1/4*norm_eta;
      Delta = 1/4*Delta;
   elseif rho > 3/4 && (stop_inner == 2 || stop_inner == 1),
      trstr = 'TR+';
      %Delta = min(2*norm_eta,Delta_bar);
      Delta = min(2*Delta,Delta_bar);
   end

   % choose new iterate based on performance
   oldgradx = fgradx;
   if rho > rho_prime,
      accept = true;
      accstr = 'acc';
      x    = x_prop;
      fx   = fx_prop;
      fgradx = fns.fgrad(x);
      norm_grad = sqrt(fns.g(x,fgradx,fgradx));
   else
      accept = false;
      accstr = 'REJ';
   end
      
   % ** Testing for Stop Criteria
   % min_outer is the minimum number of inner iterations
   % before we can exit. this gives randomization a chance to
   % escape a saddle point.
   if norm_grad < epsilon && (~useRand || k > min_outer),
      stop_outer = 1;
   end

   % stop clock for this outer step
   this_time = toc;
   % update total time
   total_time = total_time + this_time;

   % ** Record data:
   curstat.k        = k;
   curstat.ng       = norm_grad;
   curstat.fx       = fx;
   curstat.rho      = rho;
   curstat.rhonum      = rhonum;
   curstat.rhoden      = rhoden;
   curstat.time     = this_time;
   curstat.accepted = accept;
   curstat.numinner = numit;
   curstat.Delta    = Delta;
   if fns.havedist,
      curstat.dist = fns.dist(x,xsol);
   end
   if useRand,
      curstat.cauchy = used_cauchy;
   end
   stats(end+1) = curstat;

   % ** Display:
   if verbosity == 1,
      fprintf(['%3s %3s   k: %5d     num_inner: %5d     ',...
               'f: %e   |grad|: %e   %s\n'],...
              accstr,trstr,k,numit,fx,norm_grad,srstr);
   elseif verbosity > 1,
      if useRand && used_cauchy,
         fprintf('USED CAUCHY POINT\n');
      end
      fprintf(['%3s %3s    k: %5d     num_inner: %5d     %s\n'],...
              accstr,trstr,k,numit,srstr);
      fprintf('       f(x) : %e     |grad| : %e\n',fx,norm_grad);
      fprintf('      Delta : %f          |eta| : %e\n',Delta,norm_eta);
      if small_rhonum ~= 0,
         fprintf('VERY SMALL rho_num: %e\n',small_rhonum);
      else
         fprintf('        rho : %e\n',rho);
      end
      if fns.havedist,
         fprintf('       dist : %e\n',stats(end).dist);
      end
      fprintf('       Time : %f\n',this_time);
   end
   if debug > 0,
      fprintf('DBG: cos ang(eta,gradf): %d\n',testangle);
      if rho==0
         keyboard;
      end
   end
   
   % stop after max_outer iterations
   if k >= max_outer,
      if (verbosity > 0),
         fprintf('\n*** timed out -- k == %d***\n',k);
      end
      stop_outer = 1;
   end 

end  % of TR loop (counter: k)
if (verbosity > 1) || (debug > 0),
   fprintf('************************************************************************\n');
end
if (verbosity > 0) || (debug > 0)
    fprintf('Total time is %f\n',total_time);
end

varargout{1} = x;
if nargout > 1,
   varargout{2} = stats;
end
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  get_string %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret,opts] = get_string(opts,argname,argdesc,def,options)

    % Process inputs and do error-checking 
    errstr = sprintf('%s opts.%s must be: \n',argdesc,argname);
    errstr = [errstr, sprintf('%s  ',options{:})];
    if isfield(opts,argname)
        ret = getfield(opts,argname);
        valid = 0;
        if isstr(ret),
            for i = 1:length(options),
                if isequal(ret,options{i}),
                    valid = 1;
                    break;
                end
            end
        end
        if ~valid,
            error(errstr);
        end

        % remove field from opts
        opts = rmfield(opts,argname);
    else
        ret = def;
    end
    return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  get_int %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret,opts] = get_int(opts,argname,argdesc,def,lb,ub)

    if nargin < 6
        ub = inf;
    end
    if nargin < 5,
        lb = -inf;
    end
    if nargin < 4,
        mandatory = 1;
        errstr = sprintf('%s opts.%s is mandatory, an integer in [%d,%d]',...
                         argdesc,argname,lb,ub);
    else
        mandatory = 0;
        errstr = sprintf('%s opts.%s must be an integer in [%d,%d]',...
                         argdesc,argname,lb,ub);
    end

    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = getfield(opts,argname);
        valid = 0;
        % check that it is an int
        if isnumeric(ret),
            ret = floor(ret);
            % check size (1 by 1) and bounds
            if isequal(size(ret),[1 1]) && lb <= ret && ret <= ub,
                valid = 1;
            end
        end
        if ~valid,
            error(errstr);
        end

        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory == 0,
        ret = def;
    else
        error(errorstr);
    end
    return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  get_float %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret,opts] = get_float(opts,argname,argdesc,def,lb,ub)

    if nargin < 6
        ub = inf;
    end
    if nargin < 5,
        lb = -inf;
    end
    if nargin < 4,
        mandatory = 1;
        errstr = sprintf('%s opts.%s is mandatory, a scalar in [%d,%d]',...
                         argdesc,argname,lb,ub);
    else 
        mandatory = 0;
        errstr = sprintf('%s opts.%s must be a scalar in [%d,%d]',...
                         argdesc,argname,lb,ub);
    end

    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = getfield(opts,argname);
        valid = 0;
        % check that it is an int
        if isnumeric(ret),
            ret = double(ret);
            if lb <= ret && ret <= ub,
               valid = 1;
            end
        end
        if ~valid,
            error(errstr);
        end

        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory == 0,
        ret = def;
    else
        error(errstr);
    end
    return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  get_iterate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret,opts] = get_iterate(opts,argname,argdesc,def);

    if nargin < 4,
        mandatory = 1;
        errstr = sprintf('%s opts.%s is mandatory',argdesc,argname);
    else 
        mandatory = 0;
    end

    % Process inputs and do error-checking 
    if isfield(opts,argname)
        ret = getfield(opts,argname);
        % remove field from opts
        opts = rmfield(opts,argname);
    elseif mandatory==0,
        ret = def;
    else
        error(errstr);
    end
    return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% truncated CG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,inner_it,stop_tCG] = tCG(fns,x,grad,eta,Delta,theta,kappa,min_inner,max_inner,useRand,debug);
% tCG - Truncated (Steihaug-Toint) Conjugate-Gradient
% minimize <eta,grad> + .5*<eta,Hess(eta)>
% subject to <eta,eta> <= Delta^2

   % all terms involving the trust-region radius will utilize an inner product
   % w.r.t. the preconditioner; this is because the iterates grow in
   % length w.r.t. the preconditioner, guaranteeing that we will not 
   % re-enter the trust-region
   % 
   % the following recurrences for Prec-based norms and inner 
   % products come from CGT2000, pg. 205, first edition
   % below, P is the preconditioner
   % 
   % <eta_k,P*delta_k> = beta_k-1 * ( <eta_k-1,P*delta_k-1> + alpha_k-1 |delta_k-1|^2_P )
   % |delta_k|^2_P = <r_k,z_k> + beta_k-1^2 |delta_k-1|^2_P
   % 
   % therefore, we need to keep track of 
   % 1)   |delta_k|^2_P 
   % 2)   <eta_k,P*delta_k> = <eta_k,delta_k>_P
   % 3)   |eta_k  |^2_P
   % 
   % initial values are given by:
   %    |delta_0|_P = <r,z>
   %    |eta_0|_P   = 0
   %    <eta_0,delta_0>_P = 0
   % because we take eta_0 = 0

   if useRand, % and therefore, fns.haveprec == 0
      % eta (presumably) ~= 0 was provided by the caller   
      %r = grad+fns.fhess(x,eta);
      r = tgtveclincomb(1, grad, 1, fns.fhess(x, eta));
      e_Pe = fns.g(x,eta,eta);
   else % and therefore, eta == 0
      % eta = 0*grad;
      r = grad;
      e_Pe = 0;
   end
   r_r = fns.g(x,r,r);
   norm_r = sqrt(r_r);
   norm_r0 = norm_r;

   % precondition the residual
   if fns.haveprec == 0,
      z = r;
   else
      z = fns.prec(x,r);
   end
   % compute z'*r
   z_r = fns.g(x,z,r);
   d_Pd = z_r;

   % initial search direction
   % delta  = -z;
   delta = tgtveclincomb(-1, z);
   
   if useRand, % and therefore, fns.haveprec == 0
      e_Pd = fns.g(x,eta,delta);
   else % and therefore, eta == 0
      e_Pd = 0;
   end

   % pre-assume termination b/c j == end
   stop_tCG = 5;

   % begin inner/tCG loop
   j = 0;
   for j = 1:max_inner,

      Hxd = fns.fhess(x,delta);

      % compute curvature
      d_Hd = fns.g(x,delta,Hxd);

      % DEBUGGING: check that <d,Hd> = <Hd,d>
      if debug > 1,
         Hd_d = fns.g(x,Hxd,delta);
         fprintf('DBG: |d_Hd - Hd_d| (abs/rel): %e/%e\n',abs(d_Hd-Hd_d),abs((d_Hd-Hd_d)/d_Hd));
      end

      alpha = z_r/d_Hd;
      % <neweta,neweta>_P = <eta,eta>_P + 2*alpha*<eta,delta>_P + alpha*alpha*<delta,delta>_P
      e_Pe_new = e_Pe + 2.0*alpha*e_Pd + alpha*alpha*d_Pd;

      if debug > 2,
         fprintf('DBG:   (r,r)  : %e\n',r_r);
         fprintf('DBG:   (d,Hd) : %e\n',d_Hd);
         fprintf('DBG:   alpha  : %e\n',alpha);
      end

      % check curvature and trust-region radius
      if d_Hd <= 0 || e_Pe_new >= Delta^2,
         % want
         %  ee = <eta,eta>_prec,x
         %  ed = <eta,delta>_prec,x
         %  dd = <delta,delta>_prec,x
         tau = (-e_Pd + sqrt(e_Pd*e_Pd + d_Pd*(Delta^2-e_Pe))) / d_Pd;
         if debug > 2,
            fprintf('DBG:     tau  : %e\n',tau);
         end
         %eta = eta + tau*delta;
         eta = tgtveclincomb(1, eta, tau, delta);
         if d_Hd <= 0,
            stop_tCG = 1;     % negative curvature
         else
            stop_tCG = 2;     % exceeded trust region
         end
         break;
      end

      % no negative curvature and eta_prop inside TR: accept it
      e_Pe = e_Pe_new;
      %eta = eta + alpha*delta;
      eta = tgtveclincomb(1, eta, alpha, delta);

      % update the residual
      % r = r + alpha*Hxd;
      r = tgtveclincomb(1, r, alpha, Hxd);
      % re-tangentalize r
      r = fns.proj(x,r);

      % compute new norm of r
      r_r = fns.g(x,r,r);
      norm_r = sqrt(r_r);

      % check kappa/theta stopping criterion
      if j >= min_inner && norm_r <= norm_r0*min(norm_r0^theta,kappa)
         % residual is small enough to quit
         if kappa < norm_r0^theta,
             stop_tCG = 3;  % linear convergence
         else
             stop_tCG = 4;  % superlinear convergence
         end
         break;
      end

      % precondition the residual
      if fns.haveprec == 0,
         z = r;
      else
         z = fns.prec(x,r);
      end

      % save the old z'*r
      zold_rold = z_r;
      % compute new z'*r
      z_r = fns.g(x,z,r);

      % compute new search direction
      beta = z_r/zold_rold;
      %delta = -z + beta*delta;
      delta = tgtveclincomb(-1, z, beta, delta);
      % update new P-norms and P-dots
      e_Pd = beta*(e_Pd + alpha*d_Pd);
      d_Pd = z_r + beta*beta*d_Pd;

   end  % of tCG loop
   inner_it = j;
   return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Debugging testing of grad/Hess/R %%%%%%%%%%%%%%%%%%%
function TestGH(x,fns)

   % test hessian for symmetry
   for i=1:20,
      T1 = fns.randT(x);
      T2 = fns.randT(x);
      t1ht2 = fns.g( x, T1, fns.fhess(x,T2));
      ht1t2 = fns.g( x, fns.fhess(x,T1), T2);
      fprintf('<T1,H[T2]>: %13e \t <H[T1],T2>: %13e\n',t1ht2,ht1t2);
   end

   fprintf('\n');

   % test that hessian domain is tangent plane
   for i=1:20,
      T = fns.randT(x);
      HT = fns.fhess(x,T);
      PHT = fns.proj(x,HT);
      HTminusPHT = tgtveclincomb(1, HT, -1, PHT);
      errnrm = fns.g(x,HTminusPHT,HTminusPHT);
      errnrm = sqrt(errnrm);
      fprintf('|H[T] - Proj H[T]: %e\n',errnrm);
   end

   fprintf('\n');

   % test model for agreement with fhat
   for i=1:5,
      eta = fns.randT(x);
      eta = tgtveclincomb(1/sqrt(fns.g(x,eta,eta)), eta);
      for t = 10.^[0:-1:-6],
         fxe = fns.f( fns.R(x,tgtveclincomb(t, eta)));
         teta = tgtveclincomb(t, eta);
         mxe = fns.f( x) + fns.g(x, teta, fns.fgrad(x)) + .5*fns.g(x, teta, fns.fhess(x, teta));
         %f(R_x(t\eta)) - ( f(x) + <gradf(x),t\eta> + 1/2 <\eta,H(x)[\eta]> )
         fprintf('t: %13e     |m-foR|: %13e\n',t,fxe-mxe);
      end
      % move to a new point
      x = fns.R(x,eta);
   end

   fprintf('\n');

   % check that retraction is first-order approximation of the exponential
   if fns.havedist,
      xe = fns.R(x,0*eta);
      fprintf('dist(x, R(x,0)): %13e\n',fns.dist(x,xe));
   end

   fprintf('\n');

   % check that retraction is second-order approximation of the exponential
   % this is not necessary, but may be worth debugging for some retractions
   if fns.havedist,
      for i=1:5,
         eta = fns.randT(x);
         eta = eta / sqrt(fns.g(x,eta,eta));
         for t = 10.^[0:-1:-6],
            xe = fns.R(x,t*eta);
            fprintf('t: %13e     dist(x, R(x,t*eta)): %13e\n',t,fns.dist(x,xe));
         end
         % move to a new point
         x = fns.R(x,eta);
      end
   end

   return;

