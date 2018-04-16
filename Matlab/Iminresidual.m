% [x,ret] = Iminresidual(K,PL,PU,x0,itn,epsilon)

% local NLP on iDGP 
% using fminunc, perform local descent from x0 on obj fun 
%   min sum_ij (max(||x_i-x_j||^2-PU_ij^2,0)+max(PL_ij^2-||x_i-x_j||^2,0))
% where x,x0 in R^K, itn is the max number of iterations of the solver,
% epsilon is used for tolerance, and ret includes solution stats

function [x,ret] = Iminresidual(K,PL,PU,x0,itn,epsilon)
  opt = optimset();
  if (nargin < 5)
    % leave default iteration limit
  else
    opt = optimset(opt, 'MaxIter', itn);
  end
  if (nargin < 6)
    % leave default tolerances
  else      
    opt = optimset(opt, 'TolFun', epsilon, 'TolX', epsilon);
  end
  opt = optimset(opt, 'GradObj', 'off');
  [x,f] = fminunc(@(y)Ifresidual(K,PL,PU,y), x0, opt);  
  mx = Imde(PL,PU,x);
  lx = Ilde(PL,PU,x);
  %fprintf('minresidual: cpu=%.2f mde=%.2f lde=%.2f obj=%.4f\n',mx,lx,f);  
  ret = struct('realization',x,'obj',f);
end