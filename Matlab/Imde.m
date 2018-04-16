% mx = Imde(PL,PU,x)
%
% compute the interval mean distance error of a realization x

function mx = Imde(PL,PU,x)
  [n,n] = size(PL);
  mx = 0;
  m = 0;
  e = 0;
  for i=1:n
    for j=1:n
      if (PL(i,j)>0)
        Dij = norm(x(:,i)-x(:,j));
        e = max(PL(i,j)-Dij,0) + max(Dij-PU(i,j),0);
        mx = mx + e;
        m = m+1;
      end
    end
  end
  mx = mx/m;
end
