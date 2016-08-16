function x = basicSolution(A, b)
% let's assume that A is a 
[Q,R] = qr(A);
x = inv(R) * Q' * b;
end