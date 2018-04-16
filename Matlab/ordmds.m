% [x,err] = ordmds(M,K,epsilon)
%
% solve the ordinal MDS problem using PCA on the given distance 
% matrix M in R^K using the error tolerance epsilon. 
% x is the embedding and err its error

function [x,err] = ordmds(M,K,epsilon)
  if nargin < 3
    epsilon = 1e-6;
  end
  [n,n] = size(M);
  x = pca(M,K,epsilon);
  err = orderr(M,x);
end