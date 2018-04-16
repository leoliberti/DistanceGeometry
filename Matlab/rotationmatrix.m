% rotation matrix taking x to ||x|| y / ||y||
%   assumption: with x\not=0 and y\not\in span(x)
%   sequence of 2D (scaled) Givens rotation matrices
%   the matrix R returned by this function does not map x to y
%   if you want R such that R*x=y, look at alignrealization.m
% http://math.stackexchange.com/questions/525276/rotation-matrix-in-arbitrary-dimension-to-align-vector

function [R] = rotationmatrix(x,y)
  [K,n] = size(x);
  w = x/norm(x);
  z = y/norm(y);
  R = eye(K);
  for j = 2 : K
    R = givens(w,z,1,j) * R;
    G = givens2(w([1,j]),z([1,j]));
    w([1,j]) = G*w([1,j]);
  end
end
