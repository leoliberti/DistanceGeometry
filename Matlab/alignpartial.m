% [z,R,t] = alignpartial(x,y,m)
%
%    z is the best alignment of x to y, namely y,z should be as well
%    superposed as possible
%    z(:,1:K) = R y(:,1:K) + t is the best alignment of x(:,1:K) to
%    y(:,1:K), then z(:,K+1:m) = R y(:,K+1:m) + t.
%    Assumption: K < m <= |x|
%    Default: if m is missing, then m=n

function [z,R,t] = alignpartial(x,y,m)
  [K,n] = size(x);
  if (nargin < 3)
    m = n;
  end
  z = zeros(3,m);
  [regParams,z,ErrorStats] = absor(x(:,1:K),y(:,1:K));
  R = regParams.R;
  t = regParams.t;
  z(:,K+1:m) = R*x(:,K+1:m) + repmat(t, 1, m-K);
end
