% Householder reflection matrix w.r.t a hyperplane orthogonal to col vect v
function [H] = hhreflmat(v)
  [K,n] = size(v);
  w = v/norm(v);
  H = eye(K) - 2*w*w';
end
