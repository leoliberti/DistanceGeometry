% e = lde(P,x)
%
% compute max distance error between realization x and the pEDM P

function e = lde(P,x)
  [K,n] = size(x);
  e = 0;
  m = 0;
  derr = 0;
  for i=2:n
    for j=1:i-1
      if P(j,i) >= 0
        derr = abs(norm(x(:,i)-x(:,j)) - P(j,i));
        if derr > e
          e = abs(norm(x(:,i)-x(:,j)) - P(j,i));
        end
        m = m + 1;
      end
    end
  end
end
