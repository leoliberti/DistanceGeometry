% e = Irlzerror(PL,PU,x)
%
% compute average distance error between realization x and the interval
% pEDM pair [PL,PU]

function e = Irlzerror(PL,PU,x)
  [K,n] = size(x);
  e = 0;
  m = 0;
  d = 0;
  for i=2:n
    for j=1:i-1
      if PL(j,i) >= 0 & PU(j,i) >= 0
        d = norm(x(:,i)-x(:,j));
        %fprintf('i=%d, j=%d, d=%f, PL=%f, PU=%f\n', i,j,d,PL(j,i),PU(j,i));
        e = e + max(0,PL(j,i)-d) + max(0,d-PU(j,i));
        m = m + 1;
      end
    end
  end
  e = e/m;
end