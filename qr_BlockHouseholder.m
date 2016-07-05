% Block Householder QR
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 250 ed. 4)
function Q = qr_BlockHouseholder(A, r)
[m,n] = size(A);
Q = eye(m);
lambda = 1;
k = 0;

while lambda <= n
    tau = min(lambda+r-1, n);
    k = k + 1;
    
    % Now we upper triangularize A(lambda:m, lambda:tau)
    % generating Householder matrices H_lambda_1, ..., H_tau
    [v_store, beta_store] = qr_Householder_basic(A(lambda:m, lambda: tau) );
    
    % Now compute the block representation.
    [W, Y] = blockRepresentation(v_store, beta_store);
    K = W * Y';
    [rows,cols] = size(K);
    A(lambda:m, tau + 1:n) = (eye(rows,cols) - K)' * A(lambda:m, (tau + 1):n); 
    Q(:, lambda:m) = Q(:, lambda:m) * (eye(rows,cols) - W * Y');
    lambda = tau + 1;   
end

end