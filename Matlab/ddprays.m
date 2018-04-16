% [x,ret] = ddprays(K,P,U,epsilon)
%
%   solve U-based DDP approximation of SDP to find a rlz of the partial 
%   EDM P in R^K; epsilon is used to go from Gram matrix to realization.
%   U parametrizes the DDP (default: I_n)
%   THIS VERSION IMPLEMENTS THE DD CONSTRAINT USING THE EXTREME RAYS
%      OF THE DD CONE
%   x = the realization
%   ret: a return data structure with the following fields:
%     realization = x
%     error = the error between P and eucldist(x) over P ~= -1
%     rank = rank of the Gram matrix corresponding to completion of P
%     eigenvalues = eigenvalues of the Gram matrix
%     diagnostic = diagnostic messages from solver  

function [x,ret] = ddprays(K,P,U,epsilon)
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
  Z = sdpvar(n,n);
  d = sdpvar(n,1);
  dp = sdpvar(n,n);
  dm = sdpvar(n,n);

  F = [Y == U'*Z*U, d(:) >= 0, dp(:) >= 0, dm(:) >= 0];
  obj = 0;
  for i = 1:n-1
    for j = i+1:n
      if (P(i,j) ~= -1)
        F = F + [Y(i,i) + Y(j,j) - 2*Y(i,j) >= P(i,j)^2];
        obj = obj + (Y(i,i) + Y(j,j) - 2*Y(i,j));
      end
    end
  end
  
  % DDP constraint on Z 
  for i = 1:n
    for j = i:n
        if i == j
            F = F + [Z(i,i) == d(i,1) + sum(dp(i,[1:i-1 i+1:end]) + dm(i,[1:i-1 i+1:end]))]; 
        else
            F = F + [Z(i,j) == dp(i,j) - dm(i,j)];
        end
    end
  end
  
  diag = optimize(F,obj,opts);
  G = double(Y);
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
      'eigenvalues', eig(G), 'diagnostic', diag, ...
      'sol', G, 'obj', ov, 'd', double(d), 'dp', double(dp), 'dm', double(dm));
end

