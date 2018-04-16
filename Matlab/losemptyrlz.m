% y = losemptyrlz(x,save);
%
% Sometimes the PDB encoding of integer indices to atoms also includes
% non-atomic entities, such as TER. We do not want to renumber the atoms,
% so this results into entire empty columns in the distance matrices, 
% which is fine but yields realizations with spurious identically 
% zero columns. This function returns a realization without the
% identically zero columns.
% WARNING: if some points in x are really at zero, encode their indices
% in save

function y = losemptyrlz(x,save)
  if nargin < 2
    save = [];
  end
  [K,n] = size(x);
  d = 0;
  j = 1;
  for i=1:n
    if sum(x(:,i)) ~= 0 | ismember(i,save)
      y(:,j) = x(:,i);
      j = j + 1;
    end
  end
end