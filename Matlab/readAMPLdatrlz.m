% x = readAMPLdatrlz(K,filename)
%
%    read a realization stored into an AMPL file _rlz.dat or _sol.dat
function x = readAMPLdatrlz(K,fn)
  assert(K>0 && K<4, 'readAMPLdatrlz expects K in {1,2,3}');
  if K == 1
    [x(:,1)] = textread(fn, '%*d %f', 'headerlines', 1, 'whitespace', ' \b\t;');
  elseif K == 2
    [x(:,1),x(:,2)] = textread(fn, '%*d %f %f', 'headerlines', 1, 'whitespace', ' \b\t;');
  elseif K == 3    
    [x(:,1),x(:,2),x(:,3)] = textread(fn, '%*d %f %f %f', 'headerlines', 1, 'whitespace', ' \b\t;');
  end
  x = x';
end
