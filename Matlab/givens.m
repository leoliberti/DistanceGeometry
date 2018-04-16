% (scaled) Givens rotation matrix mapping x([i,j]) to y([i,j]) 
function [G] = givens(x,y,i,j)
  [K,n] = size(x);
  G = eye(K);
  G2 = givens2(x([i,j]),y([i,j]));
  G(i,i) = G2(1,1);
  G(i,j) = G2(1,2);
  G(j,i) = G2(2,1);
  G(j,j) = G2(2,2);
end
