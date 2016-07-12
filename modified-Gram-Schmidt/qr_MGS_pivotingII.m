function [Q, R, P] = qr_MGS_pivotingII(A)

[m, n] = size(A);
Q = zeros(m,n);
R = zeros(n);
P = 1 : 1 : n;

for i = 1 : n
    column_norms(i) = norm(A(:,i), 2);
end

for k=1:n
    
    % Swap columns
    [~, jmax] = max(column_norms);
    
    temp = A(:, jmax);
    A(:,jmax) = A(:,k);
    A(:,k) = temp;
    
    temp = column_norms(jmax);
    column_norms(jmax) = column_norms(k);
    column_norms(k) = temp;
    
    temp = P(jmax);
    P(jmax) = P(k);
    P(k) = temp;
    
    R(k,k) = norm(A(:,k));
    Q(:,k) = A(:,k)/R(k,k);
    R(k,k+1:n) = Q(:,k)'*A(:,k+1:n);
    A(:,k+1:n) = A(:,k+1:n) - Q(:,k)*R(k,k+1:n);
    
    % Update column norms
    for i = 1 : n
        column_norms(i) = norm(A(:,i), 2);
    end
    
end