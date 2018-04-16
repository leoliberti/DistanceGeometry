% [P,K,c,cinv] = readAMPLdatpedm(filename)
%
%    read a partial distance matrix stored into an AMPL file _gph.dat
%    return the partial EDM P and the dimension K. Also return a mapping
%    c:index->vertexlabel and cinv:vertexlabel->index
function [P,K,c,cinv] = readAMPLdatpedm(fn)
  fid = fopen(fn);
  if fid == -1
      disp('ERROR: fid is -1 (do not use symlink dirs on windows!)')
      P = 0;
      K = 0;
      return
  end
  tline = 'ciao';
  while (length(strfind(tline, 'param Kdim := ')) == 0)
    tline = fgets(fid);
  end
  B = sscanf(tline, 'param Kdim := %u;');
  K = B(1);
  m = 1;
  while (length(strfind(tline,'param : E : c I :=')) == 0)
    tline = fgets(fid);
    m = m + 1;
    if m >= 1000
      disp('cmd param : E : c I := not found');
      return;
    end
  end
  %c = containers.Map('KeyType', 'uint64', 'ValueType', 'uint64');
  %cinv = containers.Map('KeyType', 'uint64', 'ValueType', 'uint64');
  c = containers.Map('KeyType', 'double', 'ValueType', 'double');
  cinv = containers.Map('KeyType', 'double', 'ValueType', 'double');
  m = 0;
  n = 1;
  while (1)
    tline = fgets(fid);
    if (length(strfind(tline, ';')) == 0)
      B = sscanf(tline, '%u %u %f %u');
      if (length(B) == 4)
        m = m+1;
        A(m,:) = B;
        u = A(m,1);
        v = A(m,2);
        if ~ isKey(c,u)
            c(u) = n;
            cinv(n) = u;
            n = n+1;
        end
        if ~ isKey(c,v)
            c(v) = n;
            cinv(n) = v;
            n = n+1;
        end
      end
    else
      break;
    end
  end  
  n = n - 1;
  P = eye(n)-ones(n,n);
  for i=1:m
    P(c(A(i,1)),c(A(i,2))) = A(i,3);
    P(c(A(i,2)),c(A(i,1))) = A(i,3);
  end
  fclose(fid);
end

  
  