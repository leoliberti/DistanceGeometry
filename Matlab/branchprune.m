% X = branchprune(K,P,tolerance,order,x0)
%
%   compute set X of all incongruent solutions in R^K of the DDGP 
%   given by the partial matrix P with vertex order 'order', starting 
%   from the partial embedding x0; which, if not given explicitly, 
%   is automatically computed using realizeclique.m
%
%   WARNING: assumes graph has DDGP or DMDGP order

%% Function definition
function X = branchprune(K,P,tolerance,order,x0)
  t = cputime;
  
  %% user configuration
  % print debug messages
  debug = 0;      
  % print graphviz output for the BP search tree
  graphviz = 1;   
  % symmBP: stop at the first solution, 
  % obtain all others via partial reflections
  % IMPORTANT: don't activate on non-DMDGP instances
  symm = 1;       
  % node choice strategy:
  %   1 = breadth first
  %   2 = depth first
  %   3 = random
  %   4 = most promising bound first
  nodechoice = 2;

  if nargin < 3
    % no tolerance is given, use 1e-3
    tolerance = 0.001;
  end 
  %% Initialization
  X = [];
  bptol = tolerance;
  [m,n] = size(P);
  assert(m == n, 'partial matrix P must be square');
  if nargin < 4
    % no order nor starting point given, default order is:
    order = 1:n;
  else 
    % a vertex order is explictly given
    [m_0,n_0] = size(order);
    assert(n_0 == n, 'order should have same number of elements as P has columns');
    % reorder P
    Q = P;
    P = Q(order,order);
  end
  assert(isdmdgp(K,P) == 1, 'instance does not have a dmdgp order'); 
  if nargin < 5
    % starting point not given, compute it
    x0 = [realizeclique(P(1:K,1:K));zeros(1,K)];
  end
  [K, x0n] = size(x0);
  assert(x0n < n, 'x0 must have fewer columns than size(P)');
  
  % partition P into discretization and pruning distances
  [Ds,Pr] = dscrprnmat(P,K);

  % symmBP: compute the pruning group
  if symm == 1
    GP = pruninggroup(K,P);
  end
  
  %% initialize tree data structure
  level = x0n;
  bptree = tree(x0);
  rootnode = 1;
  leaves = [rootnode];
  solutions = 0;
  % pruned nodes counter
  prunednode = 1;
  
  outf = 1;
  if graphviz == 1
    outf = fopen('~/bptree.dot', 'w+t');
    if outf < 0
      outf = 1;
    end
    fprintf(outf,'// BP tree (GraphViz) by branchprune.m: K=%d, n=%d\n', K,n);
    fprintf(outf,'strict digraph bptree {\n');
    for k=1:K-1
      fprintf(outf, ' %d [color=gray];\n', k);
      fprintf(outf, ' %d -> %d [color=gray];\n', k,k+1);
    end
    fprintf(outf,' %d;\n', rootnode+K-1);
  end

  %% Main loop
  while size(leaves) > 0

    %% Node choice strategy
    if nodechoice == 1
      node = leaves(1); % earliest node still in queue - breadth first
    elseif nodechoice == 2
      node = leaves(size(leaves,2)); % latest node in queue - depth first on minus side
    elseif nodechoice == 3
      node = leaves(ceil(rand(1)*size(leaves,2))); % random leaf
    elseif nodechoice == 4
      %sdp-based (branch-and-bound)? TODO
    end
    
    %% Get current node data
    pathtoroot = fliplr(bptree.findpath(node,rootnode));
    level = size(pathtoroot,2) + x0n;
    x = cell2mat(arrayfun(@bptree.get, pathtoroot, 'uniformoutput', false));
    leaves = leaves(leaves ~= node);
    
    %% Compute position(s) for next point
    [xp,xm] = bpnext(P,Ds,Pr,x,level,bptol,debug);
    % different solutions?
    twosol = 1;
    if norm(xp-xm) < bptol
      % just one solution, essentially prune the - subtree
      xm = xp;
      twosol = 0;
    end
    
    % test xp,xm against pruning distances
    % x+
    if (xp(1) ~= Inf)
      okflag = 1;  
      for i=1:level-1
        prdst = Pr(i,level);
        if (prdst >= 0 && abs(norm(x(:,i)-xp) - prdst) > bptol)
          okflag = 0;
          if debug == 1
            fprintf(2,'bp: pruning node %d at level %d: ||x_i-x+||=%f, dist=%f\n', node, level, norm(x(:,i)-xp), prdst)
          end
          break;
        end
      end
      if (okflag == 1)
        [bptree, newnode] = bptree.addnode(node,xp);
        if level == n
          % one solution found
          solutions = solutions + 1;
          X(solutions,:,:) = [x,xp];
          if graphviz == 1
            fprintf(outf,' %d [color=green];\n', newnode+K-1);
            fprintf(outf,' %d -> %d;\n', node+K-1, newnode+K-1);
          end
          if symm == 1
            % symmBP, one solution found, exit main loop
            break;
          end
        else
          % not the last level, add node
          leaves = [leaves, newnode];
          if graphviz == 1
            fprintf(outf,' %d -> %d;\n', node+K-1, newnode+K-1);
          end
        end
      else
        if graphviz == 1
          fprintf(outf,' p%d [label="%f",color=red];\n', prunednode+K-1, abs(norm(x(:,i)-xp) - prdst));
          fprintf(outf,' %d -> p%d;\n', node+K-1, prunednode+K-1);
          prunednode = prunednode+1;
        end
      end
    end
    % x-
    if (twosol == 1 && xm(1) ~= Inf)
      okflag = 1;
      for i=1:level-1
        prdst = Pr(i,level);
        if (prdst >= 0 && abs(norm(x(:,i)-xm) - prdst) > bptol)
          okflag = 0;
          if debug == 1
            fprintf(2,'bp: pruning node %d at level %d: ||x_i-x-||=%f, dist=%f\n', node, level, norm(x(:,i)-xm), prdst)
          end
          break;
        end
      end
      if (okflag == 1)
        [bptree, newnode] = bptree.addnode(node,xm);
        if level == n
          % one solution found
          solutions = solutions + 1;
          X(solutions,:,:) = [x,xm];
          if graphviz == 1
            fprintf(outf,' %d [color=green];\n', newnode+K-1);
            fprintf(outf,' %d -> %d;\n', node+K-1, newnode+K-1);
          end
          if symm == 1
            % symmBP, one solution found, exit main loop
            break;
          end
        else 
          % not the last level, add nodes
          leaves = [leaves, newnode];
          if graphviz == 1
            fprintf(outf,' %d -> %d;\n', node+K-1, newnode+K-1);
          end
        end
        if graphviz == 1
          fprintf(outf,' %d -> %d;\n', node+K-1, newnode+K-1);
        end
      else
        if graphviz == 1
          fprintf(outf,' p%d [label="%f",color=red];\n', prunednode+K-1, abs(norm(x(:,i)-xm) - prdst));
          fprintf(outf,' %d -> p%d;\n', node+K-1, prunednode+K-1);
          prunednode = prunednode+1;
        end
      end
    end
  end
  
  %% post-processing
  if symm == 1
    % symmBP, compute all other solutions from X(1) and GP
    while size(GP) > 0
      v = GP(1);
      sX = size(X,1);
      c = sX;
      for i=1:sX
        c = c + 1;
        X(c,:,:) = partialreflection(squeeze(X(i,:,:)),v);
      end
      GP = GP(GP ~= v);
    end
  end
  
  %% end
  if debug == 1 && size(X) == 0
    fprintf(2,'branchprune: X=[], last solution lost at level %d\n', level);
  end
  if graphviz == 1
    fprintf(outf,'}\n');
    fclose(outf);
  end
  e = cputime - t;
  fprintf('branchprune: CPU = %d\n', e);
end