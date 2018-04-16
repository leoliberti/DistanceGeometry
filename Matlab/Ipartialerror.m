% e = Ipartialerror(PL,PU,x,y)
%
% compute distance error between realization x and y over the 
% nonnegative entries of the interval pEDM [PL,PU], and consider
% as distance error of dx and dy the sum of the errors of the 
% discrepancies of dx,dy with the corresponding entries in PL,PU,
% and, if nonzero, add |dx-dy|

function e = Ipartialerror(PL,PU,x,y)
  [K,n] = size(x);
  e = 0;
  m = 0;
  d = 0;
  Ex = 0;
  Ey = 0;
  Ed = 0;
  epsilon = 0;
  for i=2:n
    for j=1:i-1
      if PL(j,i) >= 0 & PU(j,i) >= 0
        dx = norm(x(:,i)-x(:,j));
        dy = norm(y(:,i)-y(:,j));
        Ex = max(0,PL(j,i)-dx) + max(0,dx-PU(j,i));
        Ey = max(0,PL(j,i)-dy) + max(0,dy-PU(j,i));
        Ed = abs(dx-dy);
        epsilon = (dx+dy)*1e-6;
        if (Ex + Ey > epsilon)
          e = e + Ex + Ey + Ed;
        end
        m = m + 1;
      end
    end
  end
  e = e/m;
end