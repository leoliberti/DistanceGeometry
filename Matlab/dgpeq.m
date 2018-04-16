% D = dgpeq(P,x)
function D = dgpeq(P,x)
  D = P;
  [n,n] = size(P);
  for i=1:n-1
    for j=i+1:n;
      if (P(i,j) >= 0)  
        D(i,j) = norm(x(:,[i])-x(:,[j]));
      end
    end
  end
end