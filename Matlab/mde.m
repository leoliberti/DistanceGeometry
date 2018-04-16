% e = mde(P,x)
%
% compute average distance error between realization x and the pEDM P

function e = mde(P,x)
  [K,n] = size(x);
  e = 0;
  m = 0;
  for i=2:n
    for j=1:i-1
      if P(j,i) >= 0
        e = e + abs(norm(x(:,i)-x(:,j)) - P(j,i));
        m = m + 1;
      end
    end
  end
  e = e/m;
end
