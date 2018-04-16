% x = jll(x0,epsilon)
%   Johnson-Lindenstrauss lemma projection

function x = jll(x0,epsilon,K)
  [H,n] = size(x0);
  if (nargin < 3)
    K = 4/(epsilon^2)*log(n);
  end
  Phi = (1/sqrt(K))*randn(K,H);
  x = Phi*x0;
end