% Block Householder QR with Recursion
% Copyright (c) 2016 by Pranay Seshadri
% Golub and Van Loan (page 251 ed. 4)
% Assumes that A has full column rank and nb is a positive blocking
% parameter.
function Q = qr_BlockHouseholder_Recursive(A, n, nb)

if n <= nb
    % Compute a thin QR factorization!
   [v_store, betav_store] = qr_Householder_basic(A);
end


end.gitignore