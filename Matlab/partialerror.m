% e = partialerror(P,x,y)
%
% compute distance error between realization x and y over the 
% nonnegative entries of the pEDM P

function e = partialerror(P,x,y)
  [K,n] = size(x);
  e = 0;
  m = 0;
  for i=2:n
    for j=1:i-1
      if P(i,j) >= 0
        e = e + abs(norm(x(:,i)-x(:,j)) - norm(y(:,i)-y(:,j)));
        m = m + 1;
      end
    end
  end
  e = e/m;
end