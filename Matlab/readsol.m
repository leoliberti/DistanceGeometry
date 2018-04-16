% x = readsol(K,filename)
%
%    read a realization stored into a DGSOL .sol formatted file
function x = readsol(K,fn)
  assert(K>0 && K<4, 'readAMPLdatrlz expects K in {1,2,3}');
  if K == 1
    [x(:,1)] = textread(fn, '%f', 'headerlines', 1, 'whitespace', ' \b\t;', 'commentstyle', 'shell');
  elseif K == 2
    [x(:,1),x(:,2)] = textread(fn, '%f %f', 'headerlines', 1, 'whitespace', ' \b\t;', 'commentstyle', 'shell');
  elseif K == 3    
    [x(:,1),x(:,2),x(:,3)] = textread(fn, '%f %f %f', 'headerlines', 1, 'whitespace', ' \b\t;', 'commentstyle', 'shell');
  end
  x = x';
end
