function varargout = irtr(fns,params)
% IRTR   Implicit Riemannian Trust-Region (with tCG inner solve)
%
% This performs a retraction-based, trust-region minimization of an objective function on a 
% Riemannian manifold.
%
% x = irtr(fns,params) returns the minimizer (a point on the manifold)
% [x,stats] = irtr(fns,params) also returns some statistics from the algorithm
%
% The struct fns contains the function handles necessary to perform the optimization:
%  fns.R(x,eta) : retract tangent vector eta at T_x M to the manifold
%  fns.g(x,eta,zeta) : Riemannian metric at T_x M 
%  fns.proj(x,eta) : Project eta from nearby T_x M to T_x M (used in iteration to combat round-off)
%  fns.f(x) : Compute the objective function at x
%  fns.fgrad(x) : Compute the gradient of f at x, returning a tangent vector
%  fns.fhess(x,eta) : Apply the Hessian of f, hess_x: T_x M -> T_x M
%  fns.rho(x,eta) : Compute the rho-ratio: rho(x,eta)
% Optionally, three other methods may be specified:
%  fns.rhosearch(x,eta,zeta,rhoprime) : Returns a value tau such that rho(x,eta+tau*zeta) >= rhoprime
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
%   params.x0        - An initial iterate. Because IRTR does not know what the manifold is, this is mandatory.
%   params.xsol      - For testing distance from solution using fns.dist
%   params.verbosity - Level of verbosity: 0 is silent, 1 has one-line info, 2 has detailed info
%   params.debug     - Debugging tests (0 is silent, 1 has some, 2 has a lot, 3 is more than you want to know)
%   params.max_outer - Maximum number of outer iterations (default: 100)
%   params.min_inner - Minimum number of inner iterations (default: 0). Only in effect if using randomized initial tangent vectors.
%   params.max_inner - Maximum number of inner iterations (default: inf). Recommended: dimension of manifold.
%   params.epsilon   - Outer Convergence tolerance (absolute)
%   params.kappa     - Inner kappa convergence tolerance
%   params.theta     - Inner theta convergence tolerance
%   params.rho_prime - Accept/reject ratio
%   params.testgh    - perform some simple numerical testing of gradient and Hessian
%
% The stats output is a struct array, with fields:
%  k  - the outer iteration number for the stats
%  ng - the norm of the gradient, sqrt(g(x,gradfx,gradfx))
%  fx - the current value under the objective function
%  rho - the performance ratio for the iterate
%  time - the wallclock time for the outer iteration
%  numinner - the number of inner iterations used to compute the next iterate
%  dist     - the distance from the solution
%
% See also rtr

% About: RTR - Riemannian Trust-Region
% (C) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science

% Modification history
% Version 0.1 - CGB - Tue Sep 25 2007
%   Started from rtr.m in GenRTR (Baker)

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
   error('Invalid arguments to irtr. See ''help irtr''');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check functions handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if   ~(isfield(fns,'R')     && isa(fcnchk(fns.R)    ,'function_handle')) ...
  || ~(isfield(fns,'g')     && isa(fcnchk(fns.g)    ,'function_handle')) ...
  || ~(isfield(fns,'proj')  && isa(fcnchk(fns.proj) ,'function_handle')) ...
  || ~(isfield(fns,'f')     && isa(fcnchk(fns.f)    ,'function_handle')) ...
  || ~(isfield(fns,'fgrad') && isa(fcnchk(fns.fgrad),'function_handle')) ...
  || ~(isfield(fns,'rho')   && isa(fcnchk(fns.rho)  ,'function_handle')) ...
  || ~(isfield(fns,'fhess') && isa(fcnchk(fns.fhess),'function_handle')),
  error('Invalid arguments to irtr. See ''help irtr''');
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
if isfield(fns,'rhosearch') && isa(fcnchk(fns.rhosearch),'function_handle'),
    fns.haverhosearch = 1;
else
    fns.haverhosearch = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = 'tCG';
verbosity =   1;
debug     =   0;
min_inner =   0;
max_inner = inf;
max_outer = 100;
epsilon   =   1e-6;
kappa     =   0.1;
theta     =   1.0;
rho_prime =   0.1;
testgh    = 0;
gamma = 2.0/3.0;
%x0        =  user must specify

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[verbosity,params] = get_int(params,'verbosity','Level of verbosity'                ,verbosity);
[debug    ,params] = get_int(params,'debug'    ,'Debugging tests'                   ,debug);
[max_outer,params] = get_int(params,'max_outer','Maximum number of outer iterations',max_outer,0);
[min_inner,params] = get_int(params,'min_inner','Minimum number of inner iterations',min_inner,0);
[max_inner,params] = get_int(params,'max_inner','Maximum number of inner iterations',max_inner,0);
[testgh   ,params] = get_int(params,'testgh'   ,'Test gradient and Hessian'         ,testgh);
[x0       ,params] = get_iterate(params,'x0'   ,'Initial iterate');
[xsol     ,params] = get_iterate(params,'xsol' ,'Solution for testing',[]);
if isempty(xsol),
    % must have fns.dist and xsol to make fns.dist() comparisons
    fns.havedist = 0;
end
[epsilon  ,params] = get_float(params,'epsilon'  ,'Outer Convergence tolerance'      ,epsilon,0);
[kappa    ,params] = get_float(params,'kappa'    ,'Inner kappa convergence tolerance',kappa,0,1);
[theta    ,params] = get_float(params,'theta'    ,'Inner theta convergence tolerance',theta,0);
[rho_prime,params] = get_float(params,'rho_prime','Accept/reject parameter'          ,rho_prime,0,1);
[gamma    ,params] = get_float(params,'gamma'    ,'Backtracking/Fronttracking factor'    ,gamma,0,1);

if testgh && fns.haverand,
   TestGH(x0,fns);
   varargout{1} = x0;
   if nargout > 1,
      varargout{2} = struct([]);
   end
   return;
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

% ** Record data:
curstat.k        = k;
curstat.ng       = norm_grad;
curstat.fx       = fx;
curstat.rho      = inf;
curstat.time = this_time;
curstat.numinner = 0;
if fns.havedist,
    curstat.dist = fns.dist(x,xsol);
end
stats(1) = curstat;

% ** Display:
if (verbosity == 1)
   fprintf(['   %5s                %5s     ',...
            'f: %e   |grad|: %e\n'],...
           '     ','     ',fx,norm_grad);
elseif (verbosity > 1)
   fprintf('************************************************************************\n');
   fprintf(['k: %5s     num_inner: %5s     %s\n'],...
           '','','______','______','');
   fprintf('       f(x) : %e       |grad| : %e\n',fx,norm_grad);
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
   eta = 0*fgradx;

   % solve TR subproblem
   if isequal(solver,'tCG'),
      [eta,numit,stop_inner,rho] = tCG(fns,x,fgradx,eta,rho_prime,theta,kappa,gamma,min_inner,max_inner,debug);
      srstr = tcg_stop_reason{stop_inner};
   end
   norm_eta = sqrt(fns.g(x,eta,eta));
   if debug > 0,
      testangle = fns.g(x,eta,fgradx) / (norm_eta*norm_grad);
   end

   if debug > 0,
      x_prop  = fns.R(x,eta);
      fx_prop = fns.f(x_prop);
      Hxeta = fns.fhess(x,eta);
      rhonum = fx-fx_prop;
      rhoden = -fns.g(x,fgradx,eta) - 0.5*fns.g(x,Hxeta,eta);
      actrho =   rhonum  / rhoden;
      fprintf('DBG:   tCG rho : %e\n',rho);
      fprintf('DBG:   act rho : %e / %e == %e\n',rhonum,rhoden,actrho);
      fprintf('DBG:   new f(x): %e\told f(x): %e\n',fx_prop,fx);
      if fx_prop > fx*(1+sqrt(sqrt(eps))),
          warning('IRTR:ObjectiveIncrease','Increase in objective function from tCG.\nCheck that increase is not significant.');
      end
      if rho*abs(rhoden) > abs(rhonum)*(1+sqrt(sqrt(eps))),
          warning('IRTR:RhosDoNotMatch','Computed rho after tCG is less than estimate returned from tCG.\nCheck that difference is not significant.');
      end
   end

   x = fns.R(x,eta);
   fx = fns.f(x);
   fgradx = fns.fgrad(x);
   norm_grad = sqrt(fns.g(x,fgradx,fgradx));
      
   % ** Testing for Stop Criteria
   if norm_grad < epsilon,
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
   curstat.time     = this_time;
   curstat.numinner = numit;
   if fns.havedist,
      curstat.dist = fns.dist(x,xsol);
   end
   stats(end+1) = curstat;

   % ** Display:
   if verbosity == 1,
      fprintf(['k: %5d     num_inner: %5d     ',...
               'f: %e   |grad|: %e   %s\n'],...
              k,numit,fx,norm_grad,srstr);
   elseif verbosity > 1,
      fprintf(['k: %5d     num_inner: %5d     %s\n'],...
              k,numit,srstr);
      fprintf('       f(x) : %e     |grad| : %e\n',fx,norm_grad);
      fprintf('        rho : %f          |eta| : %e\n',rho,norm_eta);
      if fns.havedist,
         fprintf('       dist : %e\n',stats(end).dist);
      end
      fprintf('       Time : %f\n',this_time);
   end
   if debug > 0,
      fprintf('DBG: cos ang(eta,gradf): %d\n',testangle);
      if testangle>0,
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
function [eta,inner_it,stop_tCG,rho] = tCG(fns,x,grad,eta,rho_prime,theta,kappa,gamma,min_inner,max_inner,debug);
% tCG - Truncated (Steihaug-Toint) Conjugate-Gradient
% solve for tangent vector eta: Hess(eta) = -grad

   % eta = 0*grad;
   r = grad;
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

   % initial search direction
   delta  = -z;

   % pre-assume termination b/c j == end
   stop_tCG = 5;

   % begin inner/tCG loop
   j = 0;
   % rho_x(0) = 1
   rho = 1;
   for j = 1:max_inner,

      Hxd = fns.fhess(x,delta);

      % compute curvature
      d_Hd = fns.g(x,delta,Hxd);

      % DEBUGGING: check that <d,Hd> = <Hd,d>
      if debug > 1,
         Hd_d = fns.g(x,Hxd,delta);
         fprintf('DBG: |d_Hd - Hd_d| (abs/rel): %e/%e\n',abs(d_Hd-Hd_d),abs((d_Hd-Hd_d)/d_Hd));
      end

      % compute step size and rho
      alpha = z_r/d_Hd;
      oldrho = rho;
      rho = fns.rho(x,eta+alpha*delta);

      if debug > 1,
         fprintf('DBG:   (r,r)  : %e\n',r_r);
         fprintf('DBG:   (d,Hd) : %e\n',d_Hd);
         fprintf('DBG:   alpha  : %e\n',alpha);
         fprintf('DBG:     rho  : %e\n',rho);
      end

      % check curvature and trust-region
      if d_Hd <= 0 || rho < rho_prime,
         % search for tau such that rho(x,tau*eta) >= rho_prime
         if fns.haverhosearch,
            % the user knows how to search
            tau = fns.rhosearch(x,eta,delta,rho_prime);
            rho = fns.rho(x,eta+tau*delta);
         else
            % we will have to search
            %keyboard;
            if d_Hd > 0,
               % no negative curvature, must have exited the trust-region
               %  rho(eta+    0*delta) > rho_prime > rho(eta+alpha*delta) 
               % binary search in (0,alpha)
               l = 0;
               h = alpha;
               [tau,rho] = BinaryRho(fns,rho_prime,x,eta,delta,l,h);
            else
               % negative curvature
               % then alpha < 0, so that rho(eta+alpha*delta) is meaningless 
               % start with eta+abs(alpha)*delta and start moving out
               tau_ = 0;
               rho_ = oldrho;
               tau = abs(alpha);
               rho = fns.rho(x,eta+tau*delta);
               % loop invariance:
               % rho(eta+tau_*delta) = rho_ >= rho_prime
               % after loop:
               % rho(eta+tau*delta)  = rho < rho_prime
               if debug > 2,
                   fprintf('rho(%e) == %e\n',tau_,oldrho);
                   fprintf('rho(%e) == %e\n',tau,rho);
               end
               while rho >= rho_prime,
                  tau_ = tau;
                  rho_ = rho;
                  tau = 2*tau;
                  rho = fns.rho(x,eta+tau*delta);
                  if debug > 2,
                      fprintf('rho(%e) == %e\n',tau,rho);
                  end
               end
               % this one is bad, last was good
               % do search now between [tau_,tau]
               l = tau_;
               h = tau;
               if debug > 2,
                   fprintf('running search: [%e,%e]\n',l,h);
               end
               [tau,rho] = BinaryRho(fns,rho_prime,x,eta,delta,l,h);
            end
         end
         eta = eta + tau*delta;
         if debug > 1
            fprintf('DBG:   (e,g)  : %e\n',fns.g(x,eta,grad) );
            fprintf('DBG:   (e,He) : %e\n',fns.g(x,eta,fns.fhess(x,eta)) );
         end
         if d_Hd <= 0,
            stop_tCG = 1;     % negative curvature
         else
            stop_tCG = 2;     % exceeded trust region
         end
         break;
      end

      % no negative curvature and eta_prop inside TR: accept it
      eta = eta + alpha*delta;

      % update the residual
      r = r + alpha*Hxd;
      % re-orthogonalize or re-tangentalize
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
      delta = -z + beta*delta;

   end  % of tCG loop
   inner_it = j;
   return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% binary rho search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tau,rho] = BinaryRho(fns,rho_prime,x,eta,delta,l,h)
   %  rho(eta+l*delta) > rho_prime > rho(eta+h*delta) 
   % binary search between [l,h]
   if l > h,
      % invalid input
      fprintf('ERROR: l > h\n');
      keyboard;
   end
   tau = (l+h)/2;
   rho = fns.rho(x,eta+tau*delta);
   % rho_prime is in (0,1)
   % reduce distance to rho_prime to 10% of rho_prime
   i = 0;
   while (rho-rho_prime) > rho_prime/10 || rho < rho_prime,
      i = i + 1;
      if i > 20, 
         % we have reduced range by 2^20 ~= 1e6
         % something is probably wrong
         keyboard;
      end
      if rho > rho_prime,
         l = tau;
      else
         h = tau;
      end
      tau = (l+h)/2;
      rho = fns.rho(x,eta+tau*delta);
   end
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
      fprintf('<T1,H[T2]>: %13e      H[T1],T2: %13e\n',t1ht2,ht1t2);
   end

   fprintf('\n');

   % test that hessian domain is tangent plane
   for i=1:20,
      T = fns.randT(x);
      HT = fns.fhess(x,T);
      PHT = fns.proj(x,HT);
      errnrm = fns.g(x,HT-PHT,HT-PHT);
      errnrm = sqrt(errnrm);
      fprintf('|H[T] - Proj H[T]: %e\n',errnrm);
   end

   fprintf('\n');

   % test model for agreement with fhat
   for i=1:5,
      eta = fns.randT(x);
      eta = eta / sqrt(fns.g(x,eta,eta));
      for t = 10.^[0:-1:-6],
         fxe = fns.f( fns.R(x,t*eta));
         mxe = fns.f( x) + fns.g(x,t*eta,fns.fgrad(x)) + .5*fns.g(x,t*eta,fns.fhess(x,t*eta));
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

