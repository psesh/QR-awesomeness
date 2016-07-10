% Block Householder QR
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 250 ed. 4)
% Coded July 5th 2016
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

function [v_store, betav_store] = qr_Householder_basic(A)

% Size of A
[m,n] = size(A);

for j = 1 : n - 1
    [v,betav] = house(A(j:m,j));
    H = (eye(m-j+1) - betav * (v * v') ); % I'd prefer not to use H{j}!
    A(j:m,j:n) =  H * A(j:m,j:n);
    if j < m
        A(j+1:m,j) = v(2:m - j + 1);
    end
    
    % Store
    v_store(:,j) = [zeros(j-1,1) ; v];
    betav_store(j) = betav;
end

end
