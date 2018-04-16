% [x,ret] = localsolve(K,P,x0,itn,epsilon)

% perform local descent from x0 on obj fun 
%   min sum_ij (||x_i-x_j||^2 - P_ij^2)
% where x,x0 in R^K, itn is the max number of iterations of the solver,
% epsilon is used for tolerance, and ret includes solution stats

function [x,ret] = localsolve(K,P,x0,solver,itn,epsilon)
  if (nargin < 4)
    solver = 'snopt';
  end
  if (nargin < 5)
    itn = 0; % unlimited
  end
  if (nargin < 6)
    epsilon = 0.0001;
  end
  opts = sdpsettings('solver', solver);

  tic;
  [n,n] = size(P);
  X = sdpvar(n,K);
  
  obj = 0;
  for i = 1:n-1
    for j = i+1:n
      if (P(i,j) ~= -1)
        term = 0;
        for h = 1:K
          term = term + (X(i,h) - X(j,h))^2;
        end
        term = term - P(i,j)^2;
        obj = obj + term^2;
      end
    end
  end

  diag = optimize([],obj,opts);
  x = double(X);
  mx = mde(P,x);
  lx = lde(P,x);
  ov = double(obj);
  cpu = toc;
  fprintf('localsolve: cpu=%.2f mde=%.2f lde=%.2f obj=%.4f\n',cpu,mx,lx,ov);  
  ret = struct('realization', x, 'diagnostic', diag,'obj', ov);
  
end