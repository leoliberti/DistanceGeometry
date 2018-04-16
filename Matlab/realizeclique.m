% x = realizeclique(P,epsilon)
% 
% Computes the realization of a clique in R^{n-1} from a nxn complete EDM

function x = realizeclique(P,epsilon)
  if nargin < 2
    epsilon = 1e-6;
  end
  [n,n] = size(P);
  % realize 2-cliques in constant time
  x = zeros(n-1,n);    % x_1 = 0
  x(1,2) = P(1,2); % x_21 = P_12
  if n > 2
    % realize 3-cliques in constant time
    x(1,3) = (P(1,2) + (P(1,3)^2-P(2,3)^2) / P(1,2))/2;
    x(2,3) = sqrt(P(1,3)^2 - x(1,3)^2);
    if n > 3
      % realize general n-cliques with n>3, based on x1=0
      A = 2*x(:,2)';
      b = norm(x(:,2))^2-P(2,n)^2+P(1,n)^2;
      for K=3:n-1
        A = [A ; 2*x(:,K)'];
        b = [b ; norm(x(:,K))^2-P(K,n)^2+P(1,n)^2];
        B = A(:,1:K-1);
        Binv = inv(B);
        aK = A(1:K-1,K);
        if all(aK<epsilon)
          % aK is a zero column, simple case
          y1 = Binv*b;
          yK = sqrt(P(1,n)^2 - norm(y1)^2);
        else
          % aK nonzero, compute discriminant (this should never
          % happen unless we start with a nontrivial congruence of
          % x1,x2,x3)
          lambda = 1;
          mu = 0;
          nu = -P(1,n)^2;
          for h=1:K-1
            for j=1:K-1
              lambda = lambda + Binv(h,j)^2 * aK(j)^2;
              mu = mu + Binv(h,j)^2 * b(h) * aK(j);
              nu = nu + Binv(h,h)^2 * b(h)^2;
            end
          end        
          discr2 = mu^2 - lambda*nu;
          discr = real(sqrt(discr2));
          %// for cliques, we just take the + solution
          yK = (mu + discr)/lambda;
          y1 = Binv*(b-yK*aK);
        end
        x(1:K,K+1) = [y1 ; yK];
      end
    end
  end
end