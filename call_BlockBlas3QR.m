% Call to Block QR column pivoting
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [Ablock, column_norms, perm] = call_BlockBlas3QR(m, n, idealnb, A)

% Initialize vectors perm and colnorms and set j = 1
j = 1;
perm = 1 : 1 : n;

% Compute column norms:
for j = 1 : n
    column_norms(j) = norm(A(:,j),2);
end

while j <= n
    nb = min(idealnb, n-j+1);
    [Ablock, column_norms, perm] = qr_BlockBlas3QR(m, n-j+1, j, nb, A(:, j:n), perm, column_norms);
    j = j + nb;
end
end