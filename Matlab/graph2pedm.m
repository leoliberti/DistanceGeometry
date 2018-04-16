% P = graph2pedm(G)
%
% Transforms a n x 3 array graph G, with weight({G(1),G(2)})=G(3),
%   into a partial Euclidean distance matrix P

function P = graph2edm(G)
  [m,t] = size(G);
  n = max(max(G(:,1)),max(G(:,2)));
  P = zeros(n);
  for i=1:m
    P(G(i,1),G(i,2)) = G(i,3);
    P(G(i,2),G(i,1)) = G(i,3);
  end
end