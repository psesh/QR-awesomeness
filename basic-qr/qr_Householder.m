% Code for computing QR factorization using Householder vectors
% Copyright (c) 2016 by Pranay Seshadri
function [Q,R] = qr_Householder(A, varargin)

if isempty(varargin)
    thin = 0;
else
    thin = 1;
end

% Size of A
[m,n] = size(A);

for j = 1 : n
    [v,betav] = house(A(j:m,j));
    % Problem is the step below! Requires me to store the submatrix
    % A(j:m, j:m)!!!
    A(j:m,j:n) = (eye(m-j+1) - betav * (v * v') ) * A(j:m,j:n);
    if j < m
        A(j+1:m,j) = v(2:m - j + 1);
    end
end
R = triu(A); % R is the upper triangular matrix of the "new" A

% Computation of Q using backward accumulation
k = min(m,n);
Q = eye(m,m);
for j = k : -1 : 1
    v = [1; A((j+1):m,j)];
    betav = 2/(1 + norm(A((j+1):m,j), 2)^2); % We get the beta's from the stored Householder vectors!
    Q(j:m,j:m) = Q(j:m,j:m) - (betav * (v*v') * Q(j:m,j:m));
end

% Making it thin!?
if thin == 1
    Q = Q(1:m, 1:n);
    R = R(1:n, 1:n);
end

end
