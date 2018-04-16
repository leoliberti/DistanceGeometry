% [Q,order,lost,invorder] = losemptyedm(P);
%
% Sometimes the PDB encoding of integer indices to atoms also includes
% non-atomic entities, such as TER. We do not want to renumber the atoms,
% so this results into entire empty columns in the distance matrices, 
% which is fine but yields realizations with spurious identically 
% zero columns. This function returns a distance matrix without the
% rows/columns which are totally negative (aside from the diagonal 0).
% The index array 'lost' contains the indices of the removed columns,
% its complement save=setdiff(1:n,lost) can be used in losemptyrlz.m
% The index arrays order and invorder encode the new order in terms
% of the old, and vice-versa: with respect to the new index variable k
% and the old index variable i, we have: order(k)=i and 
% invorder(i)=k or invorder(i)=0 if i is in lost

function [Q,order,lost,invorder] = losemptyedm(P)
  [n,n] = size(P);
  j = 1;
  k = 1;
  lost = [];
  order = [];
  invorder = [];
  for i=1:n
    row = P(i,:);
    r1 = row([1:i-1,i+1:n]);
    if all(r1<0) 
      lost(j) = i;
      j = j + 1;
      invorder(i) = 0;
    else
      order(k) = i;
      invorder(i) = k;
      k = k+1;
    end
  end
  keep = setdiff(1:n,lost);
  Q = P(keep,keep);
end