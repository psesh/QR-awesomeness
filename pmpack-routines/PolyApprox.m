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
m = 10;  n = 20;

Astar = A(:,1:m); % First select number of basis terms
[~,~,pivots_m] = qr(Astar', 'vector');  % QR column pivoting
pivots_m = pivots_m(1:m);

% So what happens if I add 1 more point?
k = 15;
Ahat = A(:, (m+1):m+k);
[~,~,pivots_k] = qr(Ahat', 'vector');
pivots_k = pivots_k(1:k);
[~,~,pivots_f] = qr(A(:, 1:(m+k) )', 'vector');  % QR column pivoting

% Figure out which pivots are common!
old_pivots = pivots_m;
new_pivots = pivots_k;
counter = 1;
for i = 1 : length(old_pivots)
    [boolean, location] = amIinarray(old_pivots(i), new_pivots);
    if(strcmp(boolean, 'true'))
        common_value(counter) = old_pivots(i);
        common_location(counter) = location;
        counter = counter + 1;
    end
end

% Form a new "A" with only the old pivots and the common ones! Let's assume
% that none of the new pivots are in the old one! --> WORST CASE
A_lsqr = A(old_pivots, 1:m+k);
C_lsqr = A(new_pivots, 1:m+k);
b = y(old_pivots);
d = y(new_pivots);
x_update = solve_constrained_LS(A_lsqr, C_lsqr, b, d);

% Now we proceed to solve the linearly constrained least squares problem!
x_old = A(old_pivots, 1:m) \ y(old_pivots); % Old estimate
x_naive = A([old_pivots, new_pivots], 1:m+k) \ y([old_pivots, new_pivots]);
x_noadaptive = A(pivots_f, 1:m+k) \ y(pivots_f);
error_update = norm(x(1:(m+k)) - x_update, 2)
error_old = norm(x(1:m) - x_old, 2)
error_new_naive = norm(x(1:(m+k)) - x_naive, 2)
error_no_adaptive = norm(x(1:(m+k)) - x_noadaptive, 2)
