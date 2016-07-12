% Call to Block QR column pivoting
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [A_large, perm] = call_BlockBlas3QR(m, n, idealnb, A_large)

% Initialize vectors perm and colnorms and set j = 1
perm = 1 : 1 : n;
tau = zeros(n,1);

j = 1;
while j <= n
    nb = min(idealnb, n-j+1);
    [A_step, perm_step, tau_step] = qr_BlockBlas3QR_pivoting(m, n-j+1, j, nb, A_large(:, j:n), perm(j:n));
    A_large(:,j:n) = A_step;
    tau(j:j+nb-1) = tau_step;
    perm(j:n) = perm_step;
    j = j + nb;
end


end
