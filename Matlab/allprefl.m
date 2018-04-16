% X = allprefl(x,GP)
%
% generate all the partial reflections of the Kxn realization x 
% (including the identity) derived from the pruning group GP

function X = allprefl(x,GP)
  [K,n] = size(x);
  X(1,:,:) = x;
  while size(GP) > 0
    v = GP(1);
    sX = size(X,1);
    c = sX;
    for i=1:sX
      c = c + 1;
      X(c,:,:) = partialreflection(squeeze(X(i,:,:)),v);
    end
    GP = GP(GP ~= v);
  end
end