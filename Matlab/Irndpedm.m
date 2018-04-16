% [x,P] = Irndpedm(K,n,U,s)
%   generate a random realization in [-U,U] in R^K and a 
%   partial nxn interval EDM [PL,PU] with given density s in [0,1]
%   and |PU-PL| < 2*eps

function [x,PL,PU] = rndpedm(K,n,U,s,eps)
  x = 2*U*(rand(K,n) + (-0.5)*ones(K,n));
  P = eucldist(x);
  PL = random('unif', max(P-eps*ones(n,n),0.0001), P, n, n);
  PU = random('unif', P, P+eps*ones(n,n), n, n);
  for i = 1 : n
     PL(i,i) = 0;
     PU(i,i) = 0;
  end
  for i = 1 : n
    for j = i+1 : n
      if (rand(1) > s) 
        P(i,j) = -1;
        P(j,i) = -1;
        PL(i,j) = -1;
        PL(j,i) = -1;
        PU(i,j) = -1;
        PU(j,i) = -1;
      end
    end
  end
end
