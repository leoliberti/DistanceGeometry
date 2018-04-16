% b = chirality(x)
%
%   compute the chirality of the K x n realization x, i.e. for v>K
%   let (a,1) be the hyperplane containing x(:,v-K:v-1), and b(v)=1 
%   if a'*x(:,v)>1, b(v)=-1 if a'x(:,v)<1, and b(v)=0 otherwise
%

function b = chirality(x,epsilon)
  if (nargin < 2)
    epsilon = 0.01;
  end
  [K,n] = size(x);
  b = zeros(1,n);
  for v = K+1:n
    a = x(:,v-K:v-1)' \ ones(K,1);
    rhs = a'*x(:,v);
    if rhs > 1 + epsilon
      b(v) = 1;
    elseif rhs < 1 - epsilon  
      b(v) = -1;
    else
      b(v) = 0;
    end
  end
end
