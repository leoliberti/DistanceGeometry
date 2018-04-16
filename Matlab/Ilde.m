% lx = Ilde(PL,PU,x)
%
% interval largest distance error of a realization x

function lx = Ilde(PL,PU,x)
  [n,n] = size(PL);
  lx = 0;
  e = 0;
  for i=1:n
    for j=1:n
      if (PL(i,j)>0)
        Dij = norm(x(:,i)-x(:,j));
        e = max(PL(i,j)-Dij,0) + max(Dij-PU(i,j),0);
        if e > lx
          lx = e;
        end
      end
    end
  end
end
