% [z,err,avgerr] = preflapprox(x,y)
%
%   matches the first K vertices, then for each v>K picks the 
%   partial reflection of x(v) closest to y(v). Reports the
%   norm error in err, and their average in avgerr

function [z,err,avgerr] = preflapprox(x,y)
  [K,n] = size(x);
  avgerr = 0;
  z = alignpartial(x,y);
  err = zeros(1,n);
  for i = K+1:n
    zerr = norm(y(:,i)-z(:,i));
    err(i) = zerr;
    w = partialreflection(z,i);
    werr = norm(y(:,i)-w(:,i));    
    if werr < zerr
      z = w;
      err(i) = werr;
    end
    avgerr = avgerr + err(i);
  end
  avgerr = avgerr / n;
end
