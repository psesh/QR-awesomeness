clear; close all; clc;

% Test Block-QR with pivoting
A = rand(50,30);
[m,n] = size(A);
perm = 1 : n;
ideal_nb = 5;

% Compute column norms
for j = 1 : n
   column_norms(j) = norm(A(:,j),2); 
end

% begin
j = 1;
while j <= n
    nb = min(ideal_nb, n-j+1);
    rowk = j;
    Ablock = qr_BlockBusingerGolub_pivoting(m, n - j + 1, nb, rowk, A(:, j:n) );
    j = j + nb;
end
    