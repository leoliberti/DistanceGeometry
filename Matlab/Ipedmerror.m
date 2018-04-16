% [avge,maxe] = Ipedmerror(PL,PU,D)
%
% compute the average and max edge error between a partial EDM D 
% and the bound pEDMs PL,PU. 

function [avge,maxe] = Ipedmerror(PL,PU,D)
  [n,n] = size(PL);
  avge = 0;
  maxe = 0;
  m = 0;
  e = 0;
  for i=1:n
    for j=1:n
      if (PL(i,j)>0)
        e = max(PL(i,j)-D(i,j),0) + max(D(i,j)-PU(i,j),0);
        avge = avge + e;
        if e > maxe
          maxe = e;
        end
        m = m+1;
      end
    end
  end
  avge = avge/m;
end
