% the isomap method
function [x] = isomap(P,K)
  x = pca(floydwarshall(P),K);
end
