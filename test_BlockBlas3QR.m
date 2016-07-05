clear; close all; clc;

% Setup with a random matrix
A = rand(12,11);
[m,n] = size(A);
idealnb = 4;

[Ablock, column_norms, perm] = call_BlockBlas3QR(m, n, idealnb, A);