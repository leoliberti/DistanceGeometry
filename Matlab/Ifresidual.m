% f = Ifresidual(K,PL,PU,x)
%
% returns the value of
%   min sum_ij (max(||x_i-x_j||^2-PU_ij^2,0)+max(PL_ij^2-||x_i-x_j||^2,0))
% where x is a K x n matrix representing the realization

function f = Ifresidual(K,PL,PU,x)

  [n,n] = size(PL);
  
  % function value at x
  f = 0;
  for i = 1:n-1
    for j = i+1:n
      if (PL(i,j) ~= -1)
        term = 0;
        for h = 1:K
          xijh = (x(h,i) - x(h,j));
          term = term + xijh*xijh;
        end
        termU = max(term - PU(i,j)*PU(i,j),0);
        termL = max(PL(i,j)*PL(i,j) - term,0);
        f = f + termL + termU;
      end
    end
  end
   
end
