% G = dist2gram(D)
%
% compute the Gram matrix G corresponding to a distance matrix D
% D must be a Euclidean distance matrix

function G = dist2gram(D)
  [m,n] = size(D);
  J = eye(n) - (1/n)*ones(n,n);
  G = (-1/2)* J * (D.^2) * J;
end
