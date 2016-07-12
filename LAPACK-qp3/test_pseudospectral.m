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
m = 30;  n = 30;
% Subset selection - Golub & Van Loan (page 295)
A_hat = A(:,1:m); % First select number of basis terms

% QR Column pivoting
[~,~,Pn] = callBlockBlas3QR(A_hat);  % MATLAB

% 2-nrom errors!
P2 = P(1:n);
xx = A_hat(P2',:) \ y(P2);
Pn2 = Pn(1:n);
xx2 = A_hat(Pn2,:) \ y(Pn2);
error_subset= norm(xx - x(1:m), 2)
error_qr = norm(xx2 - x(1:m), 2)

% Out of curiosity -- difference in stencil!
plot(Pn2, 'bo'); hold on; plot(P2, 'rx')


