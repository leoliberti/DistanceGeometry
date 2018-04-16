% [A, b] = gbulinsys(P,x,I,j,i)
%    geometric build-up linear system: given a partial realization x, where
%    x(I) is a submatrix of x, and (I union {i}) contained in 
%    {1,...,size(P)}, subtract ||y-x_j||^2=P_{ij}^2 from
%    ||y-x_h||^2=P_{ih}^2 for all h in I and return corresponding
%    linear system Ay=b where the h-th row of A is 2(x_h-x_j) and
%    the h-th component of b is ||x_h||^2-||x_j||^2-P_{ih}^2+P_{ij}^2
%    i is usually K+1 in the branch-and-prune, and the linsys aims to
%    determine y=x_i. j is the distance equation pivot index; in the
%    theoretical treatment we often write j=K for simplicity, i.e. we
%    subtract the last equation from the others

function [A,b] = gbulinsys(P,x,I,j,i)
  kI = size(I,2);
  K = size(x,1);
  J = I(I~=j);  % J=I\j
  A = zeros(kI-1, K);
  b = zeros(kI-1, 1);
  row = 1;
  for h=J
      A(row,:) = (x(:,h)-x(:,j))';
      b(row) = 0.5*(norm(x(:,h))^2 - norm(x(:,j))^2 - P(i,h)^2 + P(i,j)^2);
      row = row + 1;
  end
end