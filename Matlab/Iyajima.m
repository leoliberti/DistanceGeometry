% [x,ret] = Iyajima(K,PL,PU,epsilon,solver)
%
%   solve interval DGP SDP formulation with solver 'solver' to 
%   find a realization of the partial interval EDM [PL,PU] 
%   in R^K epsilon is used to go from Gram matrix to realization.
%   SDP model by Yajima
%   x = the realization
%   ret: a return data structure with the following fields:
%     realization = x
%     error = the error between P and eucldist(x) over P ~= -1
%     rank = rank of the Gram matrix corresponding to completion of P
%     eigenvalues = eigenvalues of the Gram matrix
%     diagnostic = diagnostic messages from SDP solver  

function [x,ret] = Iyajima(K,PL,PU,epsilon,solver)
  if (nargin < 4)
    epsilon = 0.0001;
    solver = 'mosek';
  end
  opts = sdpsettings( 'solver', solver);

  [n,n] = size(PL);

  Y = sdpvar(n,n);
  X = sdpvar(K,n);
  s = sdpvar(n,n);
  
  M = [ eye(K) X; X.' Y ];
  F = [ M >= 0 ];
  obj = 0;
  
  for i = 1:n-1
    for j = i+1:n
      if (PL(i,j) ~= -1)
        F = F + [Y(i,i) + Y(j,j) - 2*Y(i,j) - PL(i,j)^2 <= s(i,j)];
      end
      if (PL(i,j) ~= -1 & PU(i,j) ~= -1) 
        F = F + [2*(Y(i,i) + Y(j,j) - 2*Y(i,j)) - PL(i,j)^2 - PU(i,j)^2 <= s(i,j)];  
        obj = obj + s(i,j) - (Y(i,i) + Y(j,j) - 2*Y(i,j)) + PL(i,j)^2 + 2*Y(i,j);
        F = F + [s(i,j) >= 0];
      end
    end
  end
  
  diag = solvesdp(F,obj,opts);
  G = double(Y);
  [V,lambda] = eigs(G,K,'LM');
  x = real(eps2zero(lambda.^(1/2) * V', epsilon));
  derr = Ipedmerror(PL, PU, eucldist(x));
  grk = rank(G,epsilon);
  ret = struct('realization', x, 'error', derr, 'rank', grk, 'eigenvalues', eig(G), 'diagnostic', diag);
end

