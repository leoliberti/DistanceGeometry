% [a,failrank] = isddgp(K,P)
%
% a == 1 iff the row/column order of the pEDM P is a K-DMDGP order
% failrank contains the rank levels at which the order fails to be kDMDGP

function [a,failrank] = isddgp(K,P)
  [n,n] = size(P);
  a = 1;
  failrank = [];
  for i = 1:K-1
    for j = i+1:K
      if P(i,j) < 0
        a = 0;
        failrank = union(failrank,i+1);
        break;
      end
    end
  end
  if a == 1
    for j = K+1:n
      c = 0;
      for i = 1:j-1
        if P(i,j) >= 0
          c = c + 1;
          if c >= K
            break;
          end
        end
      end
      if c < K
        a = 0;
        failrank = union(failrank,j);
      end
    end
  end
end
