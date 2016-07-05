% Block Householder QR with Recursion
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 251 ed. 4)
% Assumes that A has full column rank and nb is a positive blocking
% parameter.
function [Q,R] = qr_BlockHouseholder_Recursive(A, n, nb)

if n <= nb
    % Compute a thin QR factorization!
    temp2 = A;
    [Q,R] = qr_Householder(temp2, 'thin');
    
else
    % Recursive calls!
    n1 = floor(n/2);
    temp2 = A(:,1:n1);
    [Q1,R11] = qr_BlockHouseholder_Recursive(temp2, n1, nb);
    R12 = Q1' * A(:,(n1+1):n);
    A(:,(n1 + 1):n) = A(:,(n1+1):n) - Q1*R12;
    temp3 = A(:,(n1 + 1):n);
    [Q2,R22] = qr_BlockHouseholder_Recursive(temp3, n-n1, nb);
    
    % Final stage.
    Q = [Q1, Q2];
    [~,cols] = size(R11);
    [rows, ~] = size(R22);
    R = [R11, R12; zeros(rows, cols), R22];
end

end