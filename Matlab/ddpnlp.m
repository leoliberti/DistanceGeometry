% [x,ret] = ddpnlp(K,P,itn,epsilon)
%
% DDP solution then PCA then local NLP solver

function [x,ret] = ddpnlp(K,P,itn,epsilon)
  grads = 0; % 0=use linesearch local NLP solver, 1=trust region
  if (nargin < 4)
    epsilon = 0.0001;
  end
  if (nargin < 3)
    itn = 5;
  end
  
  tic;
  [n,n] = size(P);
  %[x0,r0] = ddprealize(K,P);
  [x0,r0] = iterddp(K,P,itn);
  [x,ret] = minresidual(K,P,x0,grads);
  xmde = mde(P,x);
  xlde = lde(P,x);
  obj = ret.obj;
  cpu = toc;
  fprintf('ddpnlp: cpu=%.2f mde=%.2f lde=%.2f sdpobj=%.2f nlpobj=%.2f\n',...
      cpu,xmde,xlde,r0.obj,ret.obj);
end

