% Givens QR Method
% Copyright (c) 2016 by Pranay Seshadri
function [Q,R] = qr_givens(A)
[m,n] = size(A)
for j = 1 : n
    for i = m : -1 : j+1
        [c,s] = givens(A(i,j), A(i,j) );
        A([j i],j:n) = [c s ; -s c]' * A([j i],j:n) ;
    end
end
R = triu(A);

% Now compute Q!
% Computation of Q using backward accumulation
k = min(m,n);
Q = eye(m,m);
for j = k : -1 : 1
    v = [1; A(j+1:m,j)];
    betav = 2/(1 + norm(A(j+1:m,j), 2)^2); % We get the beta's from the stored Householder vectors!
    Q(j:m,j:m) = Q(j:m,j:m) - (betav * (v*v') * Q(j:m,j:m));
end
end
