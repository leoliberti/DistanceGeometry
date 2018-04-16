% w = reflection(x,y,z)
%
% computes the reflection operator mapping x to y and applies it to z
% also see reflectionoperator()
%
function w = reflection(x,y,z)
  [R,b] = reflectionoperator(x,y);
  w = R*(z-b) + b;
end