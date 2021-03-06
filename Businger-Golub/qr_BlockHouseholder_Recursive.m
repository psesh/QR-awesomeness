% Block Householder QR with Recursion
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 251 ed. 4)
% Coded July 5th 2016
% Assumes that A has full column rank and nb is a positive blocking
% parameter.
function [Q,R] = qr_BlockHouseholder_Recursive(A, n, nb)
[~,n] = size(A);
if n <= nb
    % Compute a thin QR factorization!
    [Q,R] = qr_Householder(A, 'thin');
    
else
    
    % Recursive calls!
    n1 = floor(n/2);
    TEMP = A(:,1:n1);
    [Q1,R11] = qr_BlockHouseholder_Recursive(TEMP, n1, nb);
    clear TEMP
    R12 = Q1' * A(:,(n1+1):n);
    A(:,(n1 + 1):n) = A(:,(n1+1):n) - Q1*R12;
    TEMP = A(:,(n1 + 1):n);
    [Q2,R22] = qr_BlockHouseholder_Recursive(TEMP, n - n1, nb);
    clear TEMP
    
    % Final stage.
    Q = [Q1, Q2];
    [~,cols] = size(R11);
    [rows, ~] = size(R22);
    R = [R11, R12; zeros(rows, cols), R22];
    
end

end