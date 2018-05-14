%function to calculate the max eigenvalue for a system with power method
function [eVal, eVect, iter] = Powereig(A, es)

%initial values for eigenvector and eigenvalue
eVect = ones(3, 1);
eVal = 1;

%iterate
iter = 0;

while(1)
    %calculate eigenvector and eigenvalue
    eValold = eVal;
    eVect = A * eVect;
    [eVal,index] = max(abs(eVect));
    eVect = eVect/eVal;
   
    iter = iter + 1;
    
    %check tolerance
    if eVal ~= 0
        error = norm(eVal - eValold)/norm(eValold);
        
        if error < es
            break
        end
    end
end