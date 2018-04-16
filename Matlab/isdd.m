% flag = isdd(M)
%
% return 0 iff M is diagonal dominant, or row index where DD test fails

function flag = isdd(M)
    epsilon = 0.0001;
    [n,m] = size(M);
    if n ~= m
        flag = -1;
        return
    end
    flag = 0;
    for i = 1:n
       diagonal = 0;
       for j = 1:n 
           if i ~= j
               diagonal = diagonal + abs(M(i,j));
           end
       end
       if M(i,i) < diagonal - epsilon
           flag = i
           return
    end
end