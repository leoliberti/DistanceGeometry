% [X,err] = alignment(P,rlz,sol,order,findorder,dmdgptimelimit)
%
% compare the alignment modulo partial reflections of a solution to
% the kDMDGP instance P and a bona fide solution rlz of P;
% X(1,:,:) contains the partial reflection of rlz closest to sol, 
% appropriately rotated and translated, and X(2,:,:) contains sol;
% err is the l_2 error between X(1,:,:) and sol.
% findorder = 0: never find any order using CPLEX (default but may get wrong answer)
% findorder = 1: find orders using CPLEX whenever needed (slow)
% dmdgptimelimit: max time allowed for CPLEX to run (default: 0=no limit)

function [X,err] = alignment(P,rlz,sol,order,findorder,dmdgptimelimit)
  if nargin < 4
    findorder = 0;
    order = 1:n;
    dmdgptimelimit = 0;
  end
  if nargin < 5
    findorder = 0;
    dmdgptimelimit = 0;
  end
  X = [];
  err = -1;
  % compute sizes
  [K,n] = size(rlz);
  [K1,n1] = size(sol);
  [nP,nP] = size(P);
  % check they match
  assert(n == n1, 'rlz and sol should have the same number of points');
  assert(K == K1, 'rlz and sol should have the same number of dimensions');
  % check whether there are empty rows/cols in P, rlz and sol
  [P2,ord2,lost2,invord2] = losemptyedm(P);
  [n2,n2] = size(P2);
  if nP > n2
    % we lost some empty columns/rows in P
    if n2 == n
      % reduced P is now the correct size
      P = P2;
      n = n2;
    elseif n2 < n
      % now rlz and sol have a larger size, try losing empties
      rlz = rlz(ord2);
      sol = sol(ord2);
      if n2 == n
        % we're OK now
        P = P2;
        n = n2;
      else  
        assert(n2 == n, '|rlz| != |P| even after losemptyedm() on both rlz and P');
      end
    else
      assert(n2 == n, '|rlz| != |P| even after losemptyedm()');      
    end
  end
  % compute kDMDGP order and reorder P,rlz,sol
  if findorder == 1
    a = isdmdgp(K,P);
    if a == 0
      [order,rnk] = dmdgporder(K,P,dmdgptimelimit);
    end
  end
  if size(order,2) == n
    P = P(order,order);
    rlz = rlz(:,order);
    sol = sol(:,order);
    % compute pruning group
    GP = pruninggroup(K,P);
    [z,err] = alignprefl(GP,rlz,sol);
    X(1,:,:) = z;
    X(2,:,:) = sol;
  else
    X = [];
    err = -1;
  end
end
