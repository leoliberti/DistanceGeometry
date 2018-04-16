% [yp,ym] = bpnext(P,Ds,Pr,x,i,epsilon,debug)
%
%    used by branchprune.m, called by branchprune.m
%    compute two intersection points of K spheres centered in 
%    x_{i-K},...,x_{i-1} with radii 
%    P_{i-K,i},...,P_{i-1,i}, where x is K x n
%    Original system: ||y-x_h|| = d_{ih} for all i-K <= h <= i-1 (y=x_i)
%    Choose pivot j to form 
%    Linear system: subtract j-th original eq from original system,
%      yielding K-1 x K linear system A y = b, where
%        A = (x_h-x_j)_{h\not=j}, 
%        b = 1/2 ( (||x_h||^2-||x_j||^2)-(d_ih^2-d_{jh}) )
%    Choose nonbasic k=N to form
%    Linear K-1 x K-1 subsystem B y_(\not N) + N y_N = b
%    Write as dictionary y_(\not N) = B^-1 b - B^-1 N y_N (**)
%    Choose original eq (indexed by g) to write
%        ||y||^2 - 2yx_g + ||x_g||^2 - d_ig^2 = 0,
%      replace y_(\not N) using (**), get quadratic equation in y_N
%    Solve, get two values, backsubstitute using (**) to find two 
%      vectors in R^K for the two positions of y=x_i
%    Pivot j, nonbasic N chosen to minimize the condition number of B^-1
%    Orig eq idx g chosen to maximize |y_N+ - y_N-|
%    Pr and Ds are obtained in the calling procedure using dscrprnmat

%% Function definition
function [yp,ym] = bpnext(P,Ds,Pr,x,i,epsilon,debug)
  
  %% Dimension of embedding space
  K = size(x,1);

  %% Find adjacent predecessors
  Iu = [1:i-1];
  It = Ds(i,Iu) > -1;
  I = Iu(It);
  assert(size(I,2)==K, 'partial matrix is neither DDGP nor DMDGP');
  
  %% Form the linear subsystem
  [A,b,B,N,condB,pivot,nonbasic] = bpbestlinsys(P,x,I,i,epsilon,debug);
  
  %% Form the quadratic subsystem
  Binv = B^(-1);
  [yNp,yNm] = bpbestquadratic(P,x,I,i,nonbasic,A,b,Binv,epsilon,debug);
  if (yNp == Inf)
    % no solution
    if debug == 1
      fprintf('bpnext: warning: no solutions to quadratic system\n');
    end
    yp = Inf*ones(K);
    ym = yp;
    return;
  end 
  
  %% check if it's really two solutions or just one: if one make them equal
  if abs(yNp-yNm) < epsilon
    yNp = (yNp + yNm)/2;
    yNm = yNp;
  end
  
  %% Backsubstitution to find yp, ym
  Binvb = Binv*b;
  BinvN = Binv*N;
  ypt = Binvb - BinvN*yNp;
  ymt = Binvb - BinvN*yNm;
  if (nonbasic == 1)
    yp = [yNp; ypt];
    ym = [yNm; ymt];
  elseif (nonbasic == K)
    yp = [ypt; yNp];
    ym = [ymt; yNm];
  else        
    yp = [ypt([1:nonbasic-1]); yNp; ypt([nonbasic:end])];
    ym = [ymt([1:nonbasic-1]); yNm; ymt([nonbasic:end])];
  end
end
  