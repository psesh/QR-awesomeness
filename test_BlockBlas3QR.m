clear; close all; clc;

% Setup with a random matrix
A = rand(6,4);
[m,n] = size(A);
idealnb = 2;

[Ablock, column_norms, perm] = call_BlockBlas3QR(m, n, idealnb, A);