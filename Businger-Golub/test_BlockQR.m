clear; close all; clc;

% Test Block-QR with pivoting
r = 4;
A = rand(8,6);
%Q = qr_BlockHouseholder(A, r);
nb = 4;
n = 11;
[Q, R] = qr_BlockHouseholder(A, nb);