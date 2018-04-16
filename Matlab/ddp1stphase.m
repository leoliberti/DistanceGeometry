% [x,ret] = ddp1stphase(K,P,U,epsilon)
%
%   solve U-based DDP approximation of SDP to find a rlz of the partial 
%   EDM P in R^K; epsilon is used to go from Gram matrix to realization.
%   U parametrizes the DDP (default: I_n)
%   This version solves a phase 1 problem where infeasibility is minimized
%   DDP idea by A.A. Ahmadi, application to DGP & code by L. Liberti
%   x = the realization
%   ret: a return data structure with the following fields:
%     realization = x
%     error = the error between P and eucldist(x) over P ~= -1
%     rank = rank of the Gram matrix corresponding to completion of P
%     eigenvalues = eigenvalues of the Gram matrix
%     diagnostic = diagnostic messages from solver  

function [x,ret] = ddp1stphase(K,P,U,epsilon)
  if (nargin < 4)
    epsilon = 0.0001;
  end
  yalmip('clear');
  %opts = sdpsettings('solver', 'mosek');
  opts = sdpsettings('solver', 'cplex', 'cplex.lpmethod', 4, 'cplex.barrier.crossover', -1);
  
  %tic;
  [n,n] = size(P);
  if nargin < 3
      U = eye(n);
  end  
  
  Y = sdpvar(n,n);
  X = sdpvar(K,n);
  Z = sdpvar(n+K,n+K);
  T = sdpvar(n+K,n+K);
  % slacks for equations
  sp = sdpvar(n,n);
  sm = sdpvar(n,n);

  M = [ eye(K) X; X.' Y ];
  Ubar = [ eye(K) zeros(K,n); zeros(n,K) U ];
  F = [ M == Ubar' * Z * Ubar, -T(:) <= Z(:) <= T(:) ];
  obj = 0;
  for i = 1:n-1
    F = F + [sp(i,i) >= 0, sm(i,i) >= 0];
    obj = obj + sp(i,i) + sm(i,i);
    for j = i+1:n
      F = F + [sp(i,j) >= 0, sm(i,j) >= 0];
      obj = obj + sp(i,j) + sm(i,j);
      if (P(i,j) ~= -1)
        F = F + [Y(i,i) + Y(j,j) - 2*Y(i,j) + sp(i,j)-sm(i,j) == P(i,j)^2];
        %obj = obj + (Y(i,i) + Y(j,j) - 2*Y(i,j));
      end
    end
  end
  F = F + [sp(n,n) >= 0, sm(n,n) >= 0];
  obj = obj + sp(n,n) + sm(n,n);
  
  % DDP constraint
  for i = 1:n+K
    F = F + [ Z(i,i) >= sum(T(i,[1:i-1 i+1:end]))];
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
      'sol', G, 'obj', ov, 'slackp', double(sp), 'slackm', double(sm));
end

