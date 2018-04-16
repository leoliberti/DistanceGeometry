% [order,rnk] = dmdgporder(K,P,timelimit)
%
% computes a kDMDGP order from the partial edm P;
% rnk is a mapping from the DMDGP order to the vertex set V 
% (so vtx i has rank rnk(i)), and order is the inverse function, 
% which maps a vertex to its rank (so order(i) has rank i) 
% WARNINGS: 
%  1. This function uses CPLEX, and runs it for <= timelimit seconds
%     if timelimit is not given or it is 0, then timelimit is not set
%  2. This function calls losemptyedm(P) and compresses
%     the empty rows/columns of P. In particular, if the
%     graph is disconnected, this simply establishes 
%     whether the first connected component has a kDMDGP order

function [order,rnk] = dmdgporder(K,P,timelimit)
  [n,n] = size(P);
  [Q,lost,ord,inv] = losemptyedm(P);
  [m,m] = size(Q);
  cols = m^2;
  rows = m + m + m*(K-1) + m*(m-K); 
  cplex = Cplex('dmdgporder');
  if nargin >= 3 && timelimit > 0
    cplex.Param.timelimit.Cur = timelimit;
  end
  cplex.Model.sense = 'minimize';
  cplex.Model.obj = zeros(1,cols); % pure feasibility problem
  cplex.Model.ctype = char(repmat('B',1,cols)); % binary vars
  cplex.Model.lb = zeros(cols,1);  % binary vars lower bounds
  cplex.Model.ub = ones(cols,1);   % binary vars upper bounds
  A = zeros(rows,cols); % prepare constraint matrix
  lhs = zeros(rows,1);    % prepare left hand sides
  rhs = zeros(rows,1);    % prepare right hand sides
  % unique_rank
  for i=1:m 
    for j=1:m
      A(i, (i-1)*m + j) = 1;
    end
    lhs(i) = 1;
    rhs(i) = 1;
  end
  % unique_vtx
  fprintf('ending unique_rank at row %d\n', m);
  for j=1:m
    for i=1:m
      A(m+j, (i-1)*m + j) = 1;
    end
    lhs(m+j) = 1;
    rhs(m+j) = 1;
  end
  % initial_clique
  fprintf('ending unique_vtx at row %d\n', 2*m);
  fprintf('starting initial_clique at row %d\n', 2*m + 1);
  for i=1:m
    for k=2:K
      A(2*m + (i-1)*(K-1) + k-1, (i-1)*m + k) = 1-k;
      for j=1:m
        if i ~= j && Q(i,j) ~= -1
          for h=1:k-1
            A(2*m + (i-1)*(K-1)+k-1, (j-1)*m + h) = 1;
          end  
        end
      end
      lhs(2*m + (i-1)*(K-1) + k-1) = 0;
      rhs(2*m + (i-1)*(K-1) + k-1) = inf;
    end
  end
  % kdmdgp
  fprintf('ending initial_clique at row %d\n', 2*m+m*(K-1));
  fprintf('starting kdmdgp at row %d\n', 2*m + m*(K-1) + 1);
  for i=1:m
    for k=K+1:m
      A(2*m+m*(K-1) + (i-1)*(m-K) + k-K, (i-1)*m + k) = -K;
      for j=1:m
        if i ~= j && Q(i,j) ~= -1
          for h=k-K:k-1
            A(2*m+m*(K-1) + (i-1)*(m-K) + k-K, (j-1)*m + h) = 1;
          end  
        end
      end
      lhs(2*m+m*(K-1) + (i-1)*(m-K) + k-K) = 0;
      rhs(2*m+m*(K-1) + (i-1)*(m-K) + k-K) = inf;
    end
  end
  fprintf('ending kdmdgp at row %d\n', 2*m+m*(K-1)+m*(m-K));
  cplex.Model.A = A;
  cplex.Model.lhs = lhs;
  cplex.Model.rhs = rhs;
  cplex.solve();
  fprintf('dmdgporder: CPLEX ret=%d: %s\n', cplex.Solution.status, cplex.Solution.statusstring);
  rnk = [];
  order = [];
  if cplex.Solution.status < 103
    rkf = [];
    rnk = zeros(1,n);
    order = zeros(1,m);
    for i=1:m
      for j=1:m
        if cplex.Solution.x((i-1)*m+j) == 1
          rkf(i) = j;
        end
      end
    end  
    for i=1:n
      if inv(i) > 0
        rnk(i) = rkf(inv(i));
        order(rkf(inv(i))) = i;
      else
        rnk(i) = 0;
      end
    end
  end
