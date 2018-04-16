% [x,ret] = sdpnlp(K,P,itn,epsilon)
%
% SDP solution then PCA then local NLP solver

function [x,ret] = sdpnlp(K,P,itn,epsilon)
  grads = 0; % 0=use linesearch local NLP solver, 1=trust region
  if (nargin < 4)
    epsilon = 0.0001;
  end
  if (nargin < 3)
    itn = 5;
  end
  
  tic;
  [n,n] = size(P);
  [x0,r0] = sdprealize(K,P);
  cpu0 = toc;
  [x,ret] = minresidual(K,P,x0,grads);
  xmde = mde(P,x);
  xlde = lde(P,x);
  cpu = toc;
  fprintf('sdpnlp: cpu=%.2f mde=%.2f lde=%.2f sdpobj=%.2f nlpobj=%.2f\n',...
      cpu,xmde,xlde,r0.obj,ret.obj);
  ret = struct('rlz',x, 'nlpret',ret, 'sdprlz',x0, 'sdpret',r0, 'cpu',cpu, 'cpu0',cpu0);
end

