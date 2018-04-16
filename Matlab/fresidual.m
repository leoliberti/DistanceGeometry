% f = fresidual(K,P,x)
%
% returns the value of sum_ij (||x_i-x_j||^2 - P_ij^2)^2
% where x is a K x n matrix representing the realization

function f = fresidual(K,P,x)

  [n,n] = size(P);
  
  % function value at x
  f = 0;
  for i = 1:n-1
    for j = i+1:n
      if (P(i,j) ~= -1)
        term = 0;
        for h = 1:K
          xijh = (x(h,i) - x(h,j));
          term = term + xijh*xijh;
        end
        term = term - P(i,j)*P(i,j);
        f = f + term*term;
      end
    end
  end
   
end
