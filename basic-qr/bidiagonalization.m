% Householder Bidiagonalization
%
% Returns the factorization U' * A * V, where the original matrix A is 
% repalced by a bidiagonal matrix of size m,n.
%
% See Chapter 5 of Golub & Van Loan.
% Code below is modified for both fat and tall matrices
%
% Copyright (c) 2016 by Pranay Seshadri
%
function [U, A, V] = bidiagonalization(A)

% Allocate memory!
[m,n] = size(A);
U = ones(m, m);
V = ones(n, n);

% Case 1: "TALL" MATRIX
if( m >= n)
    for j = 1 : n
        [v,b] = house(A(j:m,j));
        A(j:m, j:n) = ( eye(m-j+1) - b * (v * v') ) * A(j:m, j:n);
        if j <= n - 1
            A(j+1:m, j + 1) = v(2:m-j+1);
        end
        U = U * ( eye(m-j+1) - b * (v * v') ) 
        % Check
        U * A
        if j <= n - 2
            [v,b] = house(A(j, j+1:n)');
            A(j:m, j+1:n) = A(j:m, j+1:n) * ( eye(n-j) - b * v * v');
            A(j+1, j+2:n) = v(2:(n-j))';
            V = V * (eye(n-j) - (b * (v * v'))')
        end
    end
    
% Case 2: "FAT" MATRIX
else
    for j = 1 : m  
        [v,b] = house(A(j:m,j));
        A(j:m, j:n) = ( eye(m-j+1) - b * (v * v') ) * A(j:m, j:n);
        if j + 1 <= m
            A(j+1:m, j + 1) = v(2:m-j+1);
        end
        if j <= m - 1
            [v,b] = house(A(j, j+1:n)');
            A(j:m, j+1:n) = A(j:m, j+1:n) * ( eye(n-j) - b * v * v');
            A(j+1, j+2:n) = v(2:(n-j))';
        end
    end
    [v,b] = house(A(m, m+1:n)');
    A(m, m+1:n) = A(m, m+1:n) * ( eye(n-m) - b * v * v');
end


% Assembly of U and V matrices
end