% [x,P] = rnddmdgp(K,n,U,s)
%   generate a random realization in [-U,U] in R^K and a 
%   partial nxn euclidean distance matrix P with given density s in [0,1],
%   guaranteed to be a kDMDGP instance

function [x,P] = rnddmdgp(K,n,U,s)
  x = 2*U*(rand(K,n) + (-0.5)*ones(K,n));
  P = eucldist(x);
  for i = 1 : n
    for j = i+K+1 : n
      if (rand(1) > s)
        P(i,j) = -1;
        P(j,i) = -1;
      end
    end
  end
end

