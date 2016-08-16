%
% You will need Paul Constantine & David Gleich's pmpack routines for this!
%
% Copyright (c) 2016 by Pranay Seshadri
% University of Cambridge
clear; close all; clc;
% Setup.
s = parameter('Legendre', -1, 1); % parameter is Legendre
n_big = 10; % number of quadrature points
[p,w] = gaussian_quadrature(s,n_big); % Gauss-Legendre quadrature points
fun = @(x) exp(x(1)); % function -- you can replace this.
P = evaluate_ops(s,n_big,p); % Set up "P" and call A = P-transpose
W = diag(sqrt(w)); % diagonal matrix of sqrt(weights)
g = funceval(fun, p); % function eval'd at quadrature points
Aorig =  W' * P'; % the design matrix
[U, S, V] = svd(Aorig);
% % % % y = W' * g; % weighted function eval's
% % % % x = Aorig \ y; % "True" least squares solution!
% % % % m = 5; 
% % % % A = Aorig(:,1:m);
% % % % [~,~,pivots_m] = qr(A', 'vector');  % QR column pivoting
% % % % [Q,R,P] = qr(A');  % QR column pivoting
% % % % [~,S,~] = svd(Aorig); % Singular values of A
% % % % [~,Sr,~] = svd(R);
% % % % 
% % % % 
% % % % figure1 = figure;
% % % % set(gca, 'Yscale', 'log'); hold on; box on;
% % % % plot(S, 'ro', 'MarkerSize', 13, 'LineWidth', 2); hold on;
% % % % plot(Sr, 'bx', 'MarkerSize', 13, 'LineWidth', 2);
% % % % 
% % % % % Ok, lets do this:
% % % % b = y(pivots_m(1:m));
% % % % L = R';
% % % % L1 = L(1:m, :)
% % % % x = inv(L1) * Q * b

% % old_pivots = pivots_m(1:m);
% % Af = A(old_pivots, :);
% % bf = y(old_pivots);
% % x_old = Af \ bf % Old estimate
% % error_old = norm(x(1:m) - x_old, 2)
% 
% x_new = basicSolution(A, 



% % Af2 = Q' * R(:, 1:m)' ;

% % 
% % 
% % 
% % [Q, R] = qr(Af);
% % xqr = inv(R) * Q' * bf
% % error_qr = norm(x(1:m) - xqr, 2)
% % 



% [Q, R, P] = qr(A(:,1:m)');  % QR column pivoting
% x_new = inv(R(:, 1:m)) * Q' * b
% 
% x_old = A(old_pivots, 1:m) \ y(old_pivots); % Old estimate
% error_old = norm(x(1:m) - x_old, 2)



% % 
% % %% ONE LAYER ADAPTIVE!
% % k = 10;
% % [~,~,new_pivots] = qr_MGS_pivoting_custom(A(:, 1:m+k)', pivots_m);
% % new_pivots = new_pivots(1 : m+k);
% % 
% % %% BEST CASE!
% % [~, ~, ideal_pivots] = qr(A(:,1:m+k)', 'vector');
% % ideal_pivots = ideal_pivots(1:m+k);
% % 
% % %% ERROR ANALYSIS
% % x_old = A(old_pivots, 1:m) \ y(old_pivots); % Old estimate
% % x_adaptive = A(new_pivots, 1:m+k) \ y(new_pivots);
% % x_best = A(ideal_pivots, 1:m+k) \ y(ideal_pivots);
% % error_old = norm(x(1:m) - x_old, 2)
% % error_adaptive = norm(x(1:m+k) - x_adaptive, 2)
% % error_ideal = norm(x(1:m+k) - x_best, 2)
