% [x,ret] = barvinok(K,P,itn,epsilon,xsdp,rsdp)
%
%   Barvinok's "naive algorithm" based on concentration of measure
%   (optional) xsdp and rsdp: solution of the SDP (found by sdpnlp)

function [xstar,ret] = barvinok(K,P,xsdp,rsdp,itn,epsilon)
  grads = 0; % use linesearch local NLP solver
  if (nargin < 6)
    epsilon = 0.0001;
  end
  if (nargin < 5)
    itn = 5;
  end
  
  [n,n] = size(P);
  if nargin < 3
    % no SDP solution given, find it
    tic;
    [xsdp,rsdp] = sdprealize(K,P);
    %[xsdp,rsdp] = ddprealize(K,P);
    %[xsdp,rsdp] = iterddp(K,P,3);
    cpu0 = toc;
    x0mde = mde(P,xsdp);
    x0lde = lde(P,xsdp);
    G = rsdp.sol;
  else
    % SDP solution given
    cpu0 = rsdp.cpu0;
    x0mde = rsdp.sdpret.mde;
    x0lde = rsdp.sdpret.lde;
    G = rsdp.sdpret.sol;
  end
  fprintf('barvinok: itn=%d cpu=%.2f mde=%.2f lde=%.2f\n',0,cpu0,x0mde,x0lde);
  
  % Barvinok's naive algorithm
  tic;
  [V,lambda] = eig(G);  
  T = real(eps2zero(lambda.^(1/2) * V', epsilon));
  xstar = xsdp;
  xmdestar = x0mde;
  xldestar = x0lde;
  for i = 1:itn
    y = normrnd(0, 1/sqrt(K), K, n);
    x1 = y*T;
    x = minresidual(K,P,x1,grads);
    xmde = mde(P,x);
    xlde = lde(P,x);
    cpu1 = toc;
    if xmde < xmdestar
       xmdestar = xmde;
       xldestar = xlde;
       xstar = x;
       fprintf('barvinok: itn=%d cpu=%.2f mde=%.2f lde=%.2f\n',i,cpu1+cpu0,xmdestar,xldestar);
    end 
  end
  fprintf('barvinok: end cpu=%.2f mde=%.2f lde=%.2f\n',cpu1+cpu0,xmdestar,xldestar);
  ret = struct('realization',xstar,'cpu', cpu1+cpu0);
end
