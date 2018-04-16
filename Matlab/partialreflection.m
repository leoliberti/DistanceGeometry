% y = partialreflection(x,v,epsilon)
%
%     y is the partial reflection of x starting from vertex v wrt
%     the hyperplane defined by x(:,v-K),...,x(:,v-1)
%     Assumption: K < v <= |x|

function y = partialreflection(x,v,epsilon)
  if nargin < 3
    epsilon = 1e-6;
  end
  [K,n] = size(x);
  A = x(:,v-K:v-1);
  sysflag = 1;
  a = zeros(K,1);
  for i=1:K
    if norm(A(i,:) - repmat(A(i,1),1,K)) < epsilon
      % points all belong to a coordinate plane
      sysflag = 0;
      a(i) = 1;
      break;
    end
  end
  if sysflag == 1
    % this could fail if A is singular, see zero col case in realizeclique.m 
    a = A' \ ones(K,1); 
  end
  y = zeros(K,n);
  y(:,1:v-1) = x(:,1:v-1);
  for i = v:n
    y(:,i) = x(:,i) - 2*((x(:,i)'*a - 1)/(a'*a))*a;
  end
end