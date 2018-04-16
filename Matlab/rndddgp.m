% [x,P] = rndddgp(K,n,U,s)
%   generate a random realization in [-U,U] in R^K and a 
%   partial nxn euclidean distance matrix P with given density s in [0,1],
%   guaranteed to be a DDGP instance

function [x,P] = rndddgp(K,n,U,s)
  x = 2*U*(rand(K,n) + (-0.5)*ones(K,n));
  P = eucldist(x);
  for i = K+1 : n
    counter = 0; 
    for j = 1 : i-1
      if (counter + K < i-1 & rand(1) > s)
        counter = counter + 1;
        P(i,j) = -1;
        P(j,i) = -1;
      end
    end
  end
end

