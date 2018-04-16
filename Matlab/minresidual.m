% [x,ret] = minresidual(K,P,x0,grads,itn,epsilon)

% using fminunc, perform local descent from x0 on obj fun 
%   min sum_ij (||x_i-x_j||^2 - P_ij^2)^2
% where x,x0 in R^K, itn is the max number of iterations of the solver,
% epsilon is used for tolerance, and ret includes solution stats
% if grads = 1 then gradients are also required to fresidual

function [x,ret] = minresidual(K,P,x0,grads,itn,epsilon)
  opt = optimset();
  if (nargin < 4)
      grads = 0;
  end
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
  if grads == 1
    % output gradients
    opt = optimset(opt, 'GradObj', 'on');
    [x,f] = fminunc(@(y)fresidualgrad(K,P,y), x0, opt);  
  else
    % no gradients
    opt = optimset(opt, 'GradObj', 'off');
    [x,f] = fminunc(@(y)fresidual(K,P,y), x0, opt);  
  end
  mx = mde(P,x);
  lx = lde(P,x);
  %fprintf('minresidual: cpu=%.2f mde=%.2f lde=%.2f obj=%.4f\n',mx,lx,f);  
  ret = struct('realization',x,'obj',f);
end