% Test call_BlockBlas3QR
clear; close all; clc;
clear; close all; clc;
% Setup.
s = parameter('Legendre', -1, 1); % parameter is Legendre
n_big = 100; % number of quadrature points
[p,w] = gaussian_quadrature(s,n_big); % Gauss-Legendre quadrature points
fun = @(x) exp(20*x(1)); % function -- you can replace this.
P = evaluate_ops(s,n_big,p); % Set up "P" and call A = P-transpose
W = diag(sqrt(w)); % diagonal matrix of sqrt(weights)
g = funceval(fun, p); % function eval'd at quadrature points
A =  W' * P'; % the design matrix
y = W' * g; % weighted function eval's
x = A \ y; % "True" least squares solution!
[m,n] = size(A);
idealnb = 10;
[Aout, perm] = call_BlockBlas3QR(m, n, idealnb, A);
[Q,R,P] = qr(A, 'vector');

figure1 = figure; 
imagesc(log10(abs(Aout')));

figure2 = figure;
imagesc(log10(abs(A')))

figure3 = figure;
plot(P, 'bo'); hold on;
plot(perm, 'rx');

