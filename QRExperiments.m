% Testing
clear; close all; clc;
A = rand(10,5);

% Ok, but if MATLAB can compute the SVD of a fat matrix, then surely there
% is a way to bidiagonalize a fat matrix...or I just somehow look only at
% the transpose...?
[U,B,V] = bidiag(A);
[Q,R] = qr(B);
RHS = R * V';
LHS = U * Q;
[~,~,P] = qr(RHS, 'vector')
[~,~,P2] = qr(A, 'vector')
