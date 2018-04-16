% [x,PL,PU] = Irndddgp(K,n,U,s,eps)
%   generate a random realization in [-U,U] in R^K and a 
%   partial nxn euclidean interval distance matrix [PL,PU] 
%   with given density s in [0,1], guaranteed to be a iDDGP instance,
%   and such that the discretization distances are singleton intervals
%   PL,PU are randomly generated at most 2eps apart

function [x,PL,PU] = Irndddgp(K,n,U,s,eps)
  x = 2*U*(rand(K,n) + (-0.5)*ones(K,n));
  P = eucldist(x);
  PL = random('unif', max(P-eps*ones(n,n),0.0001), P, n, n);
  PU = random('unif', P, P+eps*ones(n,n), n, n);
  for i = 1 : n
     PL(i,i) = 0;
     PU(i,i) = 0;
  end
  for i = K+1 : n
    counter = 0; 
    for j = 1 : i-1
      if (counter + K < i-1 & rand(1) > s)
        counter = counter + 1;
        P(i,j) = -1;
        P(j,i) = -1;
        PL(i,j) = -1;
        PL(j,i) = -1;
        PU(i,j) = -1;
        PU(j,i) = -1;
      end
      if (counter + K >= i-1)
        PL(i,j) = P(i,j);
        PL(j,i) = P(j,i);
      end
    end
  end
end

