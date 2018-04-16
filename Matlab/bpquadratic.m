% [lambda,mu,nu] = bpquadratic(P,x,I,i,g,k,A,b,Binv)
%
% Given the system ||x_h-x_i||^2 = P(h,i)^2 to compute the next vertex 
% y=x_i given K known {x_h}'s for h in I, pick g-th equation and form
% linear system Eq(j)-Eq(g), getting Ay=b where A=2(x_h-x_g) and
% b=(||x_h||^2 - d_{h,i}^2). Pick nonbasic column k and let B be the set
% of basic columns (assume rk A = K-1). Let Binv = B^(-1). Replace
% y_B in Eq(g) and get quadratic [lambda*y_k^2 + mu*y_k + nu = 0].
% This function computes lambda, mu, nu.
function [lambda,mu,nu] = bpquadratic(P,x,I,i,g,k,A,b,Binv)
  K = size(A,2);
  Bb = Binv * b;
  N = A(:,k);
  BN = Binv * N;
  all = [1:K];
  basics = all(all~=k);
  xgB = x(basics,g);
  xgk = x(k,g);
  lambda = 1 + BN'*BN; % equal to 1+sum(BN.^2)
  mu = -2*((Bb'-xgB')*BN + xgk);
  nu = (Bb'-2*xgB')*Bb + norm(x(:,g))^2 - P(g,i)^2;
end
