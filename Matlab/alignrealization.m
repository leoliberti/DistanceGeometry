%  [z,R,t] = alignrealization(x,y)
%
%    best alignment of 2/3D realization x to y (z=Rx+t closest to y)
function [z,R,t] = alignrealization(x,y)
  [K,n] = size(x);
  [regParams,z,ErrorStats] = absor(x,y);
  R = regParams.R;
  t = regParams.t;
  for i=1:n
    zp(:,i) = reflection(z(:,1),y(:,1),z(:,i));
  end
  zq = norm(z-y);
  zpq = norm(zp-y);
  if zpq < zq
    z = zp;
  end
end

  %%% old code which doesn't work well:
  %%%   output rot/transl version of y s.t. x(1:K-1) matches y(1:K-1)
  %w = x - repmat(x(:,[1]),1,n); % translate x so that w1 is origin
  %z = y - repmat(y(:,[1]),1,n); % translate y so that z1 is origin
  %for i = 2 : K
  %  R = rotationmatrix(w(:,[i]),z(:,[i]));
  %  for j = i : n
  %    z(:,[j]) = R * z(:,[j]);
  %  end
  %end
  %z = z + repmat(x(:,[1]),1,n);
  

