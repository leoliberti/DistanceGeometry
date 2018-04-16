% [avg,stdev] = pgrpavgsize(dimensions,vertices,sparsity,samples)
%
% average and standar deviation of the pruning group 
% in function of graph sparsity

function [avg,stdev] = pgrpavgsize(K,n,sparsity,samples)
  s = zeros(1,samples);
  for i=1:samples
    [x,P]=rnddmdgp(K,n,10,sparsity); 
    s(i) = size(pruninggroup(K,P),2);
  end
  avg = mean(s);
  stdev = std(s);
end