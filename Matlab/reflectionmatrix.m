% returns the reflection matrix mapping x to ||x|| y / ||y||
%   assumptions: x,y are column vectors of length >= 2
%                x,y nonzero and y not in span(x)
%   the matrix R returned by this function does not map x to y
%   taken from http://math.stackexchange.com/questions/525276/rotation-matrix-in-arbitrary-dimension-to-align-vector
function [R] = reflectionmatrix(x,y)
  [K,m] = size(x);
  R = hhreflmat(x-y);
end
