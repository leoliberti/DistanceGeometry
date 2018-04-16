% 2x2 (scaled) Givens rotation matrix taking 2D vector x to y
function [G2] = givens2(x,y);
  w = x/norm(x);
  z = y/norm(y);
  costheta = w'*z;
  sintheta = (abs(1 - costheta^2))^(1/2);
  G2 = (norm(y)/norm(x))*[costheta, -sintheta; sintheta, costheta];
end
