% function b = countchr(str,chr)
%
% count the number of occurrences of characters in chr in string str

function b = countchr(str,chr)
  s = size(str,2);
  c = size(chr,2);
  b = 0;
  for i=1:s
    for j=1:c
      if str(i) == chr(j)
        b = b+1;
        break;
      end
    end
  end
end