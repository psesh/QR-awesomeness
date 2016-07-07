% Test call_BlockBlas3QR
clear; close all; clc;
A = rand(12,8);
[m,n] = size(A);
idealnb = 3;
[Ablock, column_norms, perm] = call_BlockBlas3QR(m, n, idealnb, A);