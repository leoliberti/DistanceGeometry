% GP = pruninggroup(K,P)
%
% find the generator indices (in K+1:n) for the pruning group
% induced by the kDMDGP pEDM P

function GP = pruninggroup(K,P)
  [n,n] = size(P);
  [Ds,Pr] = dscrprnmat(P,K);
  % pruning group generators
  nGP = [];
  for i=1:n-1
    for j=i+1:n
      if Pr(i,j) >= 0
        nGP = union(nGP, i+K+1:j);
      end
    end
  end
  GP = setdiff(K+1:n,nGP);
end