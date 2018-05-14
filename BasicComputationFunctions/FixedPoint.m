%Function that finds the roots to an equation using fixed point method
function [root,ea,iter,eaVec,rootVec]=FixedPoint(func,xr,es,maxit)

%allocate error vector
ea = 1;
eaVec = zeros(maxit,1);
rootVec = zeros(maxit,1);

%iterate
iter = 0;
while (1)
  xrold = xr;
  xr = func(xr);
  iter = iter + 1;
  rootVec(iter) = xr;
  
  if xr ~= 0
     ea = norm(xr - xrold)/norm(xrold) * 1;
     eaVec(iter) = ea;
  end
          
    
        disp(['iter= ', num2str(iter,'%02d'),...
            ', xr=(',num2str(xr','%7.5f '),...
            '), ea=', num2str(ea,'%8.3f')]);
  %check tolerance      
  if ea < es | iter >= maxit, break, end
  
end
root = xr;
