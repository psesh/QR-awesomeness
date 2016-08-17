% Comparing my version of MGS & Householder with pivoting
% with MATLAB's qr for random matrices
clear; close all; clc;
A = rand(23,60); [m,n] = size(A); % random matrix
[Q, P] = qr_MGS_pivoting(A)
%[Q2, R2, P2] = qr_Householder_pivoting(A);
[Q3, R3, P3] = qr(A,'vector');

figure1 = figure;
set(gca, 'FontSize', 14, 'LineWidth', 2); hold on; box on; grid on;
plot(P, 'bo', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MGS'); 
%plot(P2, 'ks', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'House'); 
plot(P3, 'rx', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'MATLAB-qr'); 
legend show;
xlabel('Pivot columns'); hold off;