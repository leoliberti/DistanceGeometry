% [x,ret] = Isdpnlp(K,PL,PU,itn,epsilon)
%
% SDP solution then PCA then local NLP solver

function [x,ret] = Isdpnlp(K,PL,PU,itn,epsilon)
  if (nargin < 5)
    epsilon = 0.0001;
  end
  if (nargin < 4)
    itn = 5;
  end
  
  tic;
  [n,n] = size(PL);
  [x0,r0] = Isdprealize(K,PL,PU);
  cpu0 = toc;
  [x,ret] = Iminresidual(K,PL,PU,x0);
  xmde = Imde(PL,PU,x);
  xlde = Ilde(PL,PU,x);
  cpu = toc;
  fprintf('Isdpnlp: cpu=%.2f mde=%.2f lde=%.2f sdpobj=%.2f nlpobj=%.2f\n',...
      cpu,xmde,xlde,r0.obj,ret.obj);
  ret = struct('rlz',x, 'nlpret',ret, 'sdprlz',x0, 'sdpret',r0, 'cpu',cpu, 'cpu0',cpu0);
end

