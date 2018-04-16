% [x,ret] = Ibarvinok(K,PL,PU,itn,epsilon,xsdp,rsdp)
%
% Barvinok's "naive algorithm" based on concentration of measure for iDGP
% (optional) xsdp and rsdp: solution of the SDP (found by sdpnlp)

function [xstar,ret] = Ibarvinok(K,PL,PU,xsdp,rsdp,itn,epsilon)
  % use linesearch local NLP solver
  if (nargin < 7)
    epsilon = 0.0001;
  end
  if (nargin < 6)
    itn = 5;
  end
  
  [n,n] = size(PL);
  if nargin < 4
    % no SDP solution given, find it
    tic;
    [xsdp,rsdp] = Isdprealize(K,PL,PU);
    %[xsdp,rsdp] = Iddprealize(K,PL,PU);
    %[xsdp,rsdp] = Iiterddp(K,PL,PU,3);
    cpu0 = toc;
    x0mde = Imde(PL,PU,xsdp);
    x0lde = Ilde(PL,PU,xsdp);
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
    x = Iminresidual(K,PL,PU,x1);
    xmde = Imde(PL,PU,x);
    xlde = Ilde(PL,PU,x);
    cpu1 = toc;
    if xmde < xmdestar
       xmdestar = xmde;
       xldestar = xlde;
       xstar = x;
       fprintf('Ibarvinok: itn=%d cpu=%.2f mde=%.2f lde=%.2f\n',i,cpu1+cpu0,xmdestar,xldestar);
    end 
  end
  fprintf('Ibarvinok: end cpu=%.2f mde=%.2f lde=%.2f\n',cpu1+cpu0,xmdestar,xldestar);
  ret = struct('realization',xstar,'cpu', cpu1+cpu0);
end
