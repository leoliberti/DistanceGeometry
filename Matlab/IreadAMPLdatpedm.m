% [PL,PU,K] = IreadAMPLdatpedm(filename)
%
%    read an interval partial EDM stored into an AMPL file _gph.dat
function [PL,PU,K] = IreadAMPLdatpedm(fn)
  fid = fopen(fn);
  tline = 'loopstart';
  while (length(strfind(tline, 'param Kdim := ')) == 0)
    tline = fgets(fid);
  end
  B = sscanf(tline, 'param Kdim := %u;');
  K = B(1);
  m = 1;
  while (length(strfind(tline,'param : E : c I cL cU :=')) == 0)
    tline = fgets(fid);
    m = m + 1;
    if m >= 1000
      disp('cmd "param : E : c I cL cU :=" not found');
      return;
    end
  end
  m = 0;
  n = 0;
  while (1)
    tline = fgets(fid);
    if (length(strfind(tline, ';')) == 0)
      B = sscanf(tline, '%u %u %f %u %f %f');
      if (length(B) == 6)
        m = m+1;
        A(m,:) = B;
        if (max(A(m,1),A(m,2)) > n)
          n = max(A(m,1),A(m,2));
        end
      end
    else
      break;
    end
  end  
  PL = eye(n)-ones(n,n);
  PU = eye(n)-ones(n,n);
  for i=1:m
    if A(i,5) ~= 0
      PL(A(i,1),A(i,2)) = A(i,5);
      PL(A(i,2),A(i,1)) = A(i,5);
      PU(A(i,1),A(i,2)) = A(i,6);
      PU(A(i,2),A(i,1)) = A(i,6);
    else
      PL(A(i,1),A(i,2)) = A(i,3);
      PL(A(i,2),A(i,1)) = A(i,3);
      PU(A(i,1),A(i,2)) = A(i,3);
      PU(A(i,2),A(i,1)) = A(i,3);        
    end
  end
  fclose(fid);
end

  
  