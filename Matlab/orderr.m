% err = orderr(M,x)
%
% ordinal distance error of the realization x w.r.t the distance rank
% matrix M

function err = orderr(M,x)
  D = eucldist(x);
  err = 0;
  for i=1:n-2
    for j=i+1:n-1
      if (D(i,j) < D(i,j+1) && M(i,j) > M(i,j+1)) || (D(i,j) > D(i,j+1) && M(i,j) < M(i,j+1))
        err = err + 1;
      end
    end
  end
end