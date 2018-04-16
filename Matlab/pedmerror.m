% e = pedmerror(P,D)
%
% compute the average error between the partial EDM P and the full or
% partial EDM D, limited to the components of P which are nonnegative

function e = pedmerror(P,D)
  [n,n] = size(P);
  e = 0;
  m = 0;
  for i=2:n
    for j=1:i-1
      if P(i,j) >= 0 
        e = e + abs(P(i,j) - D(i,j));
        m = m + 1;
      end
    end
  end
  e = e/m;
end
