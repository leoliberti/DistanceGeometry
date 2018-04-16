% [x,ret] = iterddp(K,P,itn,epsilon)
%
%   solve iterative DDP approximation of SDP to find a rlz of the partial 
%   EDM P in R^K; epsilon is used to go from Gram matrix to realization.
%   itn is the number of iterations (default: 1)
%   x = the realization
%   ret: a return data structure with the following fields:
%     realization = x
%     error = the error between P and eucldist(x) over P ~= -1
%     rank = rank of the Gram matrix corresponding to completion of P
%     eigenvalues = eigenvalues of the Gram matrix
%     diagnostic = diagnostic messages from solver  

function [x,ret] = iterddp(K,P,itn,epsilon)
  if (nargin < 4)
    epsilon = 0.0001;
  end
  if (nargin < 3)
    itn = 5;
  end
  
  [n,n] = size(P);
  U = eye(n);
  tic;
  for i = 1:itn
    %[x,ret] = ddprealize(K,P,U,epsilon);
    [x,ret] = ddprays(K,P,U,epsilon);
    Y = ret.sol;
    %[V,lambda] = eig(Y);
    %U = lambda.^(1/2) * V';
    %U = V';
    U = chol(Y + epsilon*eye(n));
    mx = mde(P,x);
    lx = lde(P,x);
    ov = ret.obj;
    cpu = toc;
    fprintf('iterddp: itn=%d cpu=%.2f mde=%.2f lde=%.2f obj=%.4f\n',i,cpu,mx,lx,ov);
  end
end
      