clear; close all; clc;

% Test Block-QR with pivoting
r = 4;
A = randi(12,11);
Q = qr_BlockHouseholder(A, r);