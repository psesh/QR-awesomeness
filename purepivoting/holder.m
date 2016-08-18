% Using Businger and Golub
%
% Copyright (c) 2016 by Pranay Seshadri
function perm  = qr_Householder_pivoting(A)
[m,n] = size(A);
perm = 1 : n;
column_norms = zeros(n,1); % Initialize column norms vector
Q = eye(m,m);

% Compute column norms
for j = 1 : n
    column_norms(j) = norm(A(1:m,j),2)^2;
end

% Reduction steps
for k = 1 : min(m,n) - 1
    
    % Compute max column norm
    [~,j_star] = max(column_norms(k:n));
    j_star = j_star + (k - 1);
    
    % Swap jth and pth columns!
    if(k ~= j_star)
        
        temp = A(1:m,k);
        A(1:m,k) = A(1:m,j_star);
        A(1:m,j_star) = temp;
        clear temp;
        
        temp = perm(k);
        perm(k) = perm(j_star);
        perm(j_star) = temp;
        clear temp;
        
    end
    
    % Reduction -- compute Householder matrix
    [v,betav] = house(A(k:m,k));
    H = (eye(m-k+1) - betav * (v * v') ); % I'd prefer not to use H{j}!
    A(k:m,k:n) =  H * A(k:m,k:n);
    Q(:,k:m) = Q(:,k:m) -  Q(:,k:m) * (v * v' * betav);
    
    % Update the remaining column norms
    if(k~=n)
        for j = k + 1 : n
            column_norms(j) = norm(A(1:m, j),2)^2;
        end
    end
end


R = triu(A);

end
