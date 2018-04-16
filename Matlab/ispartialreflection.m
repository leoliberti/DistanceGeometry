% [t,z,prefl,err,sumerr] = ispartialreflection(x,y)
%
%   if, for each v<=n, whenever prefl(v)=1 we apply a partial reflection
%   to x at v, and obtain y, then t=1, otherwise t=0
%   either because of a zero versus nonzero chirality (rare -- 
%   signaled by err(v) == -1), or because the reflected x(v), 
%   called z(v), is relatively too far from the given y(v) 
%   (signaled by err(v) = 1); sumerr is a cumulative continuous norm error

function [t,z,prefl,err,sumerr] = ispartialreflection(x,y,epsilon)
  if nargin < 3
    epsilon = 0.01;
  end
  [K,n] = size(x);
  meandistance = 0;
  sumerr = 0;
  for v = 2:n
    meandistance = meandistance + norm(x(:,v-1)-x(:,v));
  end
  for v = 2:n
    meandistance = meandistance + norm(y(:,v-1)-y(:,v));
  end
  meandistance = meandistance / (2*n-2);
  z = alignpartial(x,y);
  cz = chirality(z);
  cy = chirality(y);
  prefl = zeros(1,n);
  t = 1;
  err = zeros(1,n);
  for i = K+1:n
    % deal with unmatched collinearity
    if cy(i)*cz(i) == 0 && cy(i)+cz(i) ~= 0
      % not a partial reflection
      t = 0;
      err(i) = -1;
    end    
    if cy(i) == cz(i)
      % same chirality, no reflection
      prefl(i) = 0;
    else
      % chirality is different, need to reflect
      prefl(i) = 1;
      z = partialreflection(z,i);
      cz = chirality(z);      
    end
    nerr = norm(z(:,i) - y(:,i));
    sumerr = sumerr + nerr;
    if nerr > epsilon*meandistance
      % ||z_i - y_i|| > epsilon ||z_i - z_{i-1}||, reject
      t = 0;
      err(i) = 1;
    end
  end
end