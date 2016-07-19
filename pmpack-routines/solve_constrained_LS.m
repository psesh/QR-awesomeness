% Simple function that takes in a value and an array and determines
% if the value is in the array.
%
% Copyright (c) 2016 by Pranay Seshadri
% University of Cambridge
%
function x_update = solve_constrained_LS(A_lsqr, C_lsqr, b, d)
[~,R1] = qr(C_lsqr', 0); % Thin QR
[Q,~] = qr(C_lsqr');
u = inv(R1) * d;
A_cl = A_lsqr * Q;
A_cl_1 = A_cl(:, 1: length(u) );
A_cl_2 = A_cl(:, length(u) + 1 : end);
r = b - A_cl_1*u;
v = A_cl_2 \ r;
x_update = Q * [u; v]; % New estimate

end