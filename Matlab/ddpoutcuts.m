% [x,ret] = ddpoutcuts(K,P,U,epsilon)
%
%   solve U-based DDP approximation of SDP to find a rlz of the partial 
%   EDM P in R^K; epsilon is used to go from Gram matrix to realization.
%   U parametrizes the DDP (default: I_n)
%   THIS VERSION HAS THE DUAL CONE, SO IT'S AN OUTER APPROXIMATION;
%     OBSERVATION: ALWAYS GETS THE OPTIMAL VALUE AS SDP, SO CUT 
%     NEGATIVE EIGENVECTORS v USING vXv >= 0 CUTS 
%   Ideas by Tardella, Scozzari, Salgado, Liberti
%   x = the realization
%   ret: a return data structure with the following fields:
%     realization = x
%     error = the error between P and eucldist(x) over P ~= -1
%     rank = rank of the Gram matrix corresponding to completion of P
%     eigenvalues = eigenvalues of the Gram matrix
%     diagnostic = diagnostic messages from solver  

function [x,ret] = ddpraysout(K,P,U,epsilon)
  if (nargin < 4)
    epsilon = 0.0001;
  end
  yalmip('clear');
  %opts = sdpsettings('solver', 'mosek');
  %opts = sdpsettings('solver', 'cplex', 'cplex.lpmethod', 4, 'cplex.barrier.crossover', -1);
  opts = sdpsettings('solver', 'cplex');
  
  %tic;
  [n,n] = size(P);
  if nargin < 3
      U = eye(n);
  end  
  
  Y = sdpvar(n,n);
  F = [];
  obj = 0;
  for i = 1:n
      % these come from v U X U v >= 0 when v=e_i for some i and U=I
      F = F + [Y(i,i) >= 0];
  end
  for i = 1:n-1
    for j = i+1:n
      if (P(i,j) ~= -1)
        F = F + [Y(i,i) + Y(j,j) - 2*Y(i,j) >= P(i,j)^2];
        F = F + [Y(i,i) + Y(j,j) + 2*Y(i,j) >= 0];
        obj = obj + (Y(i,i) + Y(j,j) - 2*Y(i,j));
      else
        % these from vUXUv>=0 when v=e_i+-e_j for i<j and U=I
        F = F + [Y(i,i)+Y(j,j)-2*Y(i,j) >= 0, Y(i,i)+Y(j,j)+2*Y(i,j) >= 0];          
      end
    end
  end
  
  lm = -1;
  itn = 0;
  while lm < 0
    optret = optimize(F,obj,opts);
    G = double(Y);
    [V,lambda] = eig(G);
    negevs = [diag(lambda) < 0];
    Vn = V(:, diag(lambda) < 0);
    [nn,nm] = size(Vn);
    for i = 1:nm 
        F = F + [Vn(:,i)'*Y*Vn(:,i) >= 0];
    end
    itn = itn + 1;
    lm = lambda(1,1);
    fprintf('ddpoutcuts: itn=%d of=%f nr=%d evmin=%f\n',itn,double(obj),... 
            sum(negevs), lm);
  end
  
  [V,lambda] = eigs(G,K,'LM');
  x = real(eps2zero(lambda.^(1/2) * V', epsilon));
  derr = pedmerror(P, eucldist(x));
  grk = rank(G,epsilon);
  mx = mde(P,x);
  lx = lde(P,x);
  ov = double(obj);
  %cpu = toc;
  cpu = -1;
  fprintf('ddp: cpu=%.2f mde=%.2f lde=%.2f rk=%d obj=%.4f\n',cpu,mx,lx,grk,ov);  
  ret = struct('realization', x, 'error', derr, 'rank', grk, ...
      'eigenvalues', eig(G), 'diagnostic', optret, ...
      'sol', G, 'obj', ov);
end

