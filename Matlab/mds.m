% x = mds(D,epsilon)
%
%   multidimensional scaling: find realization x given 
%   complete distance matrix. The dimensionality of K is 
%   determined by the number of eigenvalues of the Gram 
%   matrix which are smaller than epsilon
function x = mds(D,epsilon)
  G = dist2gram(D);
  K = nnz(eps2zero(eig(G), epsilon));
  [V,lambda] = eigs(G,K,'LM');
  x = real(eps2zero(lambda.^(1/2) * V', epsilon));
%end function

