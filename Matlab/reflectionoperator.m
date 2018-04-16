% [R,bbar] = reflectionoperator(x,y)
%
%   returns a reflection operator mapping x to y, where
%   R is the Householder reflection matrix I - 2aa, with a=(x-y)/||x-y||,
%   and bbar=(0,...,0,b/a_j,0,...,0) with b = a'(x+y)/2 and a_j=argmin
%   component of a closest to 1.
%   The application of this operator to z is achieved by
%   reflection(x,y,z), and consist in w=R*(z-bbar)+bbar
%
%   Assumptions: x,y are column vectors of length >= 2
%                x,y nonzero and y not in span(x)

function [R,bbar] = reflectionoperator(x,y)
  [K,m] = size(x);
  a = (x-y)/norm(x-y);
  b = a'*((x+y)/2);
  R = eye(K) - 2*a*a';
  [mval,mind] = min(abs(a-1));
  bbar = zeros(K,1);
  bbar(mind) = b/a(mind);
  R = eye(K) - 2*a*a';
end
