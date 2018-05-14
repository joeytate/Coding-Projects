%Lu factorization by gaussian elimination

function [L, U] = LUfactor(A, n)

%initalize lower matrix as identity and upper as current A
L = eye(n);
U = A;

%solve upper matrix by forward elimination
for k = 1:n-1
    for i = k+1:n
        factor = U(i,k)/U(k,k);
        L(i,k) = factor; %update lower matrix with elimination factor
        U(i,:) = U(i,:) - factor * U(k,:);
    end
end

%output lower and upper matrix results
fprintf('lower matrix L:')
L

fprintf('upper matrix U:')
U