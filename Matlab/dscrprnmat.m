% [Discr,Prun] = dscrprnmat(P,K)
%    partition P into discretization and pruning distances; both 
%    are partial matrices where missing values are -1. K is the
%    dimensionality of the embedding space (i.e. =K discretization
%    distances)
%    P must be a DDGP matrix
function [Discr,Prun] = dscrprnmat(P,K)
  n = size(P,2);
  Discr = -1*ones(n,n) + eye(n);
  Prun = Discr;
  % first K-clique yields all discretization distances
  for i=1:K-1
    for j=i+1:K
      Discr(i,j) = P(i,j);
      Discr(j,i) = P(j,i);
    end
  end
  % rest of the distance matrix
  for i=K+1:n
    ddst = 0;
    for j=i-1:-1:1
      if (P(i,j)>-1)
        if (ddst < K)
          Discr(i,j) = P(i,j);
          Discr(j,i) = P(j,i);
          ddst = ddst + 1;
        else
          Prun(i,j) = P(i,j);
          Prun(j,i) = P(j,i);
        end
      end
    end
  end
end