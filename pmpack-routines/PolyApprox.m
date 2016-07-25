%
% You will need Paul Constantine & David Gleich's pmpack routines for this!
%
% Copyright (c) 2016 by Pranay Seshadri
% University of Cambridge
clear; close all; clc;
% Setup.
s = parameter('Legendre', -1, 1); % parameter is Legendre
n_big = 100; % number of quadrature points
[p,w] = gaussian_quadrature(s,n_big); % Gauss-Legendre quadrature points
fun = @(x) exp(3*x(1)); % function -- you can replace this.
P = evaluate_ops(s,n_big,p); % Set up "P" and call A = P-transpose
W = diag(sqrt(w)); % diagonal matrix of sqrt(weights)
g = funceval(fun, p); % function eval'd at quadrature points
A =  W' * P'; % the design matrix
y = W' * g; % weighted function eval's
x = A \ y; % "True" least squares solution!
m = 20;  n = 20;

[~,~,pivots_m] = qr(A(:,1:m)', 'vector');  % QR column pivoting
old_pivots = pivots_m(1:m);

%% ONE LAYER ADAPTIVE!
k = 10;
[~,~,new_pivots] = qr_MGS_pivoting_custom(A(:, 1:m+k)', pivots_m);
new_pivots = new_pivots(1 : m+k);

%% BEST CASE!
[~, ~, ideal_pivots] = qr(A(:,1:m+k)', 'vector');
ideal_pivots = ideal_pivots(1:m+k);

%% ERROR ANALYSIS
x_old = A(old_pivots, 1:m) \ y(old_pivots); % Old estimate
x_adaptive = A(new_pivots, 1:m+k) \ y(new_pivots);
x_best = A(ideal_pivots, 1:m+k) \ y(ideal_pivots);
error_old = norm(x(1:m) - x_old, 2)
error_adaptive = norm(x(1:m+k) - x_adaptive, 2)
error_ideal = norm(x(1:m+k) - x_best, 2)
