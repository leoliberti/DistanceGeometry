% [order,rank] = ddgporder(K,P)
%
% computes a DDGP order from the partial edm P;
% rank is a mapping from the DMDGP order to the vertex set V 
% (so vtx i has rank rank(i)), and order is the inverse function, 
% which maps a vertex to its rank (so order(i) has rank i) 
% WARNING: 
%   This function calls losemptyedm(P) and compresses
%   the empty rows/columns of P. In particular, if the
%   graph is disconnected, this simply establishes 
%   whether the first connected component has a kDMDGP order

function [order,rank] = ddgporder(K,P)
  order = [];
  rank = [];
  [n,n] = size(P);
  [Q,lost,ord,inv] = losemptyedm(P);
  [m,m] = size(Q);
  % cycle over all K-tuplets
  allC = nchoosek(1:m,K);
  sallC = size(allC,1);
  foundflag = 0;
  for c = 1:sallC
    cliqueflag = 1;
    for i = allC(c,:)
      for j = allC(c,:)
        if i < j && Q(i,j) < 0
          cliqueflag = 0;
          break;
        end
      end
    end
    if cliqueflag == 1
      [or2,rk2] = ddgporderfromclique(K,Q,allC(c,:));
      if size(or2) > 0
        % as soon as we find a DDGP order, get out of the loop
        foundflag = 1;
        break;
      end
    end
  end
  if foundflag == 1
    % combine DDGP order with "lose empty" order
    for i=1:n
      if inv(i) > 0
        rank(i) = rk2(inv(i));
        order(rank(i)) = i;
      else
        rank(i) = 0;
      end
    end
  end
end  
