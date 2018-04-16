% M = rndordmds(n,U)
%
% random n x n ordinal MDS instance on ranks [1,...,U]
function M = rndordmds(n,U)
  M = ceil(random('unif', 0.5001, 100*U, n,n)/100); 
  M = M - diag(diag(M)); 
  M = symmfromupper(M);
end
