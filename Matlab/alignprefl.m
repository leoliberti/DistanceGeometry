% [z,err] = alignprefl(GP,x,y)
% 
% given the pruning group GP and two Kxn realizations x,y, 
% find min_g,R,t ||g(R x + t) - y||_2 where g in GP, R is a rotation
% and t is a translation applied to x(:,1:K)

function [z,err] = alignprefl(GP,x,y)
  w = alignpartial(x,y);
  W = allprefl(w,GP);
  sW = 2^size(GP,2);
  err = +Inf;
  mi = 0;
  for i=1:sW
    t = norm(squeeze(W(i,:,:)) - y);
    if t < err
      err = t;
      mi = i;
    end
  end
  z = squeeze(W(mi,:,:));
end