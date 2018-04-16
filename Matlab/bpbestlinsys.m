% [A,b,B,N,condB,pivot,nonbasic] = bpbestlinsys(P,x,I,permI,level,epsilon,debug)
%
%    used by branchprune.m, called by bpnext.m
%    finds the linear subsystem with condition number closest to 1
%    over all pivot eqn idx j in (perm)I and nonbasic col idx in 1:K
function [A,b,B,N,condB,pivot,nonbasic] = bpbestlinsys(P,x,I,i,epsilon,debug)
  K = size(x,1);
  condB = Inf;
  for j=I  % pivot (index of equation to subtract)
    for k=1:K  % nonbasic candidate
      [At,bt,Bt,Nt,condBt] = bplinsys(P,x,I,j,i,k);
      if (abs(condBt-1) < abs(condB-1))
         A = At;
         b = bt;
         B = Bt;
         N = Nt;
         condB = condBt;
         pivot = j;
         nonbasic = k;
         %return;  % remove if you want to find j,k optimizing condition number
      end
    end        
  end
end
