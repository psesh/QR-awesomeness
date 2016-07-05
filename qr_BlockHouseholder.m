% Block Householder QR
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 250 ed. 4)
function [Q,R] = qr_BlockHouseholder()

Q = eye(m);
lambda = 1;
k = 0;

while lambda <= n
    tau = min(lambda + r - 1, n);
    k = k + 1;
    
    % Now we upper triangularize A(lambda:m, lambda:tau)
    % generating Householder matrices H_lambda_1, ..., H_tau
    HouseholderMatrices = qr_Householder_basic(A(lambda:m, lambda: tau);
    
    % Now compute the block representation.
    [W, Y] = blockRepresentation(HouseholderMatrices);
    
    
end








end