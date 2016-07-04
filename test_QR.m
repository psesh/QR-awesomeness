% Test QR algorithms
clear all; clc;

A = rand(5,5);

% MATLAB's default
[Q,R,P] = qr(A, 'vector');

% Businger and Golub approach
[R2,P2] = qr_BusingerGolub_pivoting(A);

% Modified Gram Schmidt
[Q3,R3, P3] = qr_MGS_pivoting(A)