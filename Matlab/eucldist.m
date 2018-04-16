% return the Euclidean matrix from a realization K x n matrix

function [D] = eucldist(x)
  D = sqeucldist(x).^(1/2);
end

  
