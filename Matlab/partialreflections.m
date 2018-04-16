% y = partialreflections(x,U,epsilon)
%
% apply a set of partial reflections w.r.t. vertices in U to x

function y = partialreflections(x,U,epsilon)
  if nargin < 3
    epsilon = 1e-6;
  end
  y = x;
  for v = U
    y = partialreflection(y,v,epsilon);
  end
end
