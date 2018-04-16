% C = allcliques(K,P)
%
% find all cliques of K elements in the graph described by the pEDM P

function C = allcliques(K,P)
  [n,n] = size(P);
  C = [];
  if n >= K
    Cidx = C;
    allC = nchoosek(1:n,K);
    sallC = size(allC,1);
    for c = 1:sallC
      cliqueflag = 1;
      for i = allC(c,:)
        for j = allC(c,:)
          if i < j && P(i,j) < 0
            cliqueflag = 0;
            break;
          end
        end
      end
      if cliqueflag == 1
        Cidx = union(Cidx, c);
      end
    end
    C = allC(Cidx,:);
  end 
end
