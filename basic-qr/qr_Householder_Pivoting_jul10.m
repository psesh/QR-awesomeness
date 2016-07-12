clear; close all; clc;
% Quick and dirty QR w/ column pivoting using Householder
A = [ 1 3 5 1; 2 -1 2 1; 1 4 6 1; 4 5 10 1];
Aorig = A;
[m,n] = size(A);
pivots = 1 : 1 : n; colnorm = 1 : 1 : n;

% Start with k = 1
k = 1;

% Compute the maximum column norm;
for j = 1 : n
    colnorm(j) =  norm(A(:,j), 2)^2;
end
[~, pmax] = max(colnorm);

[v, betav] = sub_function(Ak, Apmax)

H1 = eye(m,n) - betav * (v * v'); % 
A2 = H1 * A;

% Now update the remaining n-1 column norms:
for j = 2 : n
    colnorm(j) = colnorm(j) - A2(1,j).^2;
end

% Now because the norms are correctly sorted we can go straight to the
% Householder computation again. 
[v, betav] = house(A2(2:m, 2));
H2 = eye(m-1,n-1) - betav * (v * v');
A3 = H2 * A2(2:m, 2:n) ;

% Once again update all the column norms
for j = 3 : n
   colnorm(j) = colnorm(j) - A3(1,j-1).^2; 
end