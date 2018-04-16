
function [objval, lambda, X] = bisectiontest(Xin, Xout, P)

    [n,n] = size(Xin);

    eps = 0.01;
    for lambda = 0:eps:1
        X = lambda*Xin + (1-lambda)*Xout;
        d = eig(X);
        if d(1) >= 0
           break
        end
    end

    objval = 0;
    for i = 1:n-1
        for j = i+1:n
            if (P(i,j) ~= -1)
                objval = objval + (X(i,i) + X(j,j) - 2*X(i,j));
            end
        end
    end    
    
end