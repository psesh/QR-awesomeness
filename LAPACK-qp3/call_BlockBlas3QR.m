% Call to Block QR column pivoting
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [Ablock, column_norms, perm] = call_BlockBlas3QR(m, n, idealnb, A_large)

% Initialize vectors perm and colnorms and set j = 1
perm = 1 : 1 : n;

% Compute column norms:
for j = 1 : n
    column_norms(j) = norm(A(:,j),2);
end

j = 1;
while j <= n
    nb = min(idealnb, n-j+1);
    [A_step, column_norms, permutations] = qr_BlockBlas3QR(m, n-j+1, j, nb, A_large(:, j:n), perm(j:n), column_norms(j:n));
    j = j + nb;
end


end
