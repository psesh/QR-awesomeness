% test QR-MGS with pivoting
clear; close all; clc;

A = rand(4,10); [m,n] = size(A);
%[Q, R, P] = qr_MGS_pivoting(A);
[Q, R, P] = qr(A);
A_updated = A * P;
for j = 1 : n
   col_norms(j) = norm(A_updated(:,j), 2);
end
col_norms

%[Q2,R2,P2] = qr(A, 'vector');

