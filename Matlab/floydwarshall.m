% Floyd-Warshall algorithm to complete a partial Euclidean distance matrix
% assume input is square matrix
function [D] = floydwarshall(P)
  [m,n] = size(P);
  D = symmfromupper(P);
  for i = 1 : n
    for j = 1 : n
      if (D(i,j) < 0 || D(i,j) >= Inf)
        D(i,j) = Inf;
      end
    end
  end
  for i = 1 : n
    for j = 1 : n
      for k = 1 : n
        if (D(i,k) + D(k,j) < D(i,j))
          D(i,j) = D(i,k) + D(k,j);
        end
      end
    end
  end
end
