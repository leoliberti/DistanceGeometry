% A = pedm2adj(P)
%
% This returns the adjacency matrix of a partial distance matrix P
function A = pedm2inc(P)
  A = P;
  [n,n] = size(P);
  for i=1:n
    for j=1:n
      if (P(i,j) <= 0) 
        A(i,j) = 0;
      else
        A(i,j) = 1;
      end
    end
  end
end
