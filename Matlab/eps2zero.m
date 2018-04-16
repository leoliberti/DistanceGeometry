% any entry with abs value <= epsilon is turned to a zero
function [y] = eps2zero(x,epsilon)
  y = x;
  [m,n] = size(y);
  for i = 1:m
    for j = 1:n
      if (abs(y(i,j)) <= epsilon) 
        y(i,j) = 0.0;
      end
    end
  end
end
