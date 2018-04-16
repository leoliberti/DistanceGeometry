% [x,P] = rndNMRedm(K,n,U,s,w)
%   generate a random realization in [-U,U] in R^K and a 
%   partial nxn euclidean distance matrix P with given density s in [0,1];
%   then rewrite (100w)% of the given entries with wrong random values

function [x,P,E] = rndNMRedm(K,n,U,s,w)
  x = 2*U*(rand(K,n) + (-0.5)*ones(K,n));
  P = eucldist(x);
  E = zeros(n,n);
  for i = 1 : n
    for j = i+1 : n
      if (rand(1) > s) 
        P(i,j) = -1;
        P(j,i) = -1;
      elseif (rand(1) < w)
        P(i,j) = max(0,P(i,j) + P(i,j)*(rand(1)-0.5));
        P(j,i) = P(i,j);
        E(i,j) = 1;
        E(j,i) = 1;
      end
    end
  end  
end
