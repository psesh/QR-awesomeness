function [Q, R, P] = ssrqrcp(A)
[m,n] = size(A);
Q = zeros(m, m);
k = min(m,n);
R = zeros(k,n);
p = 4;
l = k;
O = randn(l,m);
B = O * A;
[Qb, Rb, P] = qr(B); size(B)
A1 = A * P;
[Q, R1] = qr( A1(:, 1:k) );
R12 = Q(:,1:k)' * A1(:, k+1:n );
end