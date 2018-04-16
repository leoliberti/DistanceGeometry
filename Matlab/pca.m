% x = pca(D,K,epsilon)
%
% principal component analysis - use K largest eigenvalues with K given
% assume D is a distance matrix; epsilon is optional. Returns embedding
% x in R^K
function x = pca(D,K,epsilon)
  if nargin < 3
    epsilon = 1e-6;
  end
  [V,lambda] = eigs(dist2gram(D),K,'LM');
  x = real(eps2zero(lambda.^(1/2) * V', epsilon));
end
