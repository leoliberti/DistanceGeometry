% [yNm, yNp] = bpbestquadratic(P,x,I,i,nonbasic,A,b,Binv,epsilon,debug)
%
%   used by branchprune.m, called by bpnext.m
%   choose the best quadratic subsystem over all eqn idx g in I
function [yNm, yNp] = bpbestquadratic(P,x,I,i,k,A,b,Binv,epsilon,debug)
  rootdiff = 0;
  for g=I 
    [lambda,mu,nu] = bpquadratic(P,x,I,i,g,k,A,b,Binv);
    mu2 = mu^2;
    omu2 = order(mu2);
    ln4 = 4*lambda*nu;
    discr2 = mu2 - ln4;
    if omu2 == order(ln4)
      epsilon = epsilon*10^omu2;
    end
    if (discr2 < -epsilon)
      % no real solution 
      yNp = Inf;
      yNm = Inf;
    else
      discr = real(sqrt(discr2));
      yNpt = (-mu + discr) / (2*lambda);
      yNmt = (-mu - discr) / (2*lambda);
      if (yNpt - yNmt >= rootdiff) 
        rootdiff = yNpt - yNmt;
        yNp = yNpt;
        yNm = yNmt;
      end
      return % as soon as we find a feasible solution, otherwise try again
    end
  end
