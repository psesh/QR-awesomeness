clear; close all; clc;

% Setup.
s = parameter('Legendre', -1, 1); % parameter is Legendre
n_big = 200; % number of quadrature points
[p,w] = gaussian_quadrature(s,n_big); % Gauss-Legendre quadrature points
fun = @(x) exp(20*x(1)); % function -- you can replace this.
P = evaluate_ops(s,n_big,p); % Set up "P" and call A = P-transpose
W = diag(sqrt(w)); % diagonal matrix of sqrt(weights)
g = funceval(fun, p); % function eval'd at quadrature points
A =  W' * P'; % the design matrix
y = W' * g; % weighted function eval's
x = A \ y; % "True" least squares solution!
m = 30;  n = 30;

Astar = A(:,1:m); % First select number of basis terms
[~,~,pivots_m] = qr(Astar', 'vector');  % QR column pivoting
pivots_m = pivots_m(1:m);

% So what happens if I add 50 more points?
k = 50;
Ahat = A(:, 1:k);
[~,~,pivots_k] = qr(Ahat', 'vector');
pivots_k = pivots_k(1:k);

counter = 1;
for i = m+1 : k
    value(counter) = cond(A(pivots_m, [1:m, i]));
    counter = counter + 1;
end

plot(value, 'bo')
