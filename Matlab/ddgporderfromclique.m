%%% LEO 150810: THIS IS WRONG - re-code using Alg 9 from introDG book
% [order,rank] = ddgporderfromclique(K,P,C)
%
% solve the TOP/DVOP/DDGPO on the instance (K,P) starting from clique C
% Assumption: C is a K-clique

function [order,rnk] = ddgporderfromclique(K,P,C)
  [n,n] = size(P);
  order = [];
  rnk = [];
  work = zeros(n-K,2);
  a = 1;
  for i=K+1:n
    c = 0;
    for j = 1:i-1;  # wrong - order never changes
      if P(j,i) > 0
        c = c - 1;
      end
    end
    work(i-K,1) = i;
    work(i-K,2) = c;
    if -c < K
      a = 0;
      break;
    end
  end
  if a == 1
    ord = sortrows(work,2)';
    order = [C ord(1,:)];
    rnk = zeros(1,n);
    for i=1:n
       rnk(order(i)) = i;
    end
  end
end
