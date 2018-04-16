% [A,b,B,N,condB] = bplinsys(P,x,I,j,i,k)
%     see gbulinsys.m for P,x,I,i,k arguments; k must be a (nonbasic) column
%     index in the first output argument A of gbulinsys() (hence an 
%     element of [1:K], where x is K times n); i is the BP level, i.e.
%     the vertex index whose position we're trying to find. 
%     This function returns
%     N = the k-th column of A (the nonbasic) and B = A with N removed.
%     Assumption: A is K-1 x K.
function [A,b,B,N,condB] = bplinsys(P,x,I,j,i,k)
  [A,b] = gbulinsys(P,x,I,j,i);
  N = A(:,k);
  B = A(:,[1:k-1,k+1:end]);
  condB = cond(B);
end