%function to find equation roots using bisection method
function [root,eaVec,iter]=Bisect(func,xl,xu,es,maxit)

%allocate vector for relative error values
ea = 1;
eaVec = zeros(maxit,1);

%iterate
xr = 0;
iter = 0;

while iter < maxit
    xrold = xr;
    xr = (xu + xl)/2;
    iter = iter + 1;
    eaVec(1) = 1;
    %find error after 2nd iteration and check tolerance
    if iter > 1
        ea = norm(xr - xrold)/norm(xrold);
        eaVec(iter) = ea;
          
        disp(['iter= ', num2str(iter,'%02d'),...
            ', xr=(',num2str(xr','%7.5f '),...
            '), ea=', num2str(ea,'%8.3f')]);
        if ea < es
            break
        end
    end
    if sign(func(xr)) == sign(func(xl))
        xl = xr;
    else
        xu = xr;
    end
end
root = xr;
