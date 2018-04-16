% s = pruningsparsity(K,P)
%
% sparsity of the pruning edges in a DDGP instance (K,P)

function s = pruningsparsity(K,P)
  [n,n] = size(P);
  [Ds,Pr] = dscrprnmat(P,K);
  nedges = 0;
  for i=1:n-1
    for j=i+1:n
      if Pr(i,j) >= 0
        nedges = nedges + 1;
      end
    end
  end
  nedges = nedges - (n - K + 1);
  s = nedges / (n*(n-1)/2.0);
end