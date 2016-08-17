% Householder Bidiagonalization
%
% Returns the factorization U' * A * V, where the original matrix A is
% repalced by a bidiagonal matrix of size m,n.
%
% See Chapter 5 of Golub & Van Loan.
% Code below is modified for both fat and tall matrices
%
% Copyright (c) 2016 by Pranay Seshadri
%
function [U, A, V] = bidiagonalization(A)

% Allocate memory!
[m,n] = size(A);
Aorig = A;
U = eye(m,m);
V = eye(n,n);

% Case 1: "TALL" MATRIX
if( m >= n)
    last_entry = n;
    for j = 1 : n
        [v,b] = house(A(j:m,j));
        
        % --- added lines ---
        U(j:m,j:m) = eye(m-j+1)  - b * v * v' * U(j:m,j:m);
        U
        U * A
        %--------------------
        A(j:m, j:n) = ( eye(m-j+1) - b * (v * v') ) * A(j:m, j:n);
        if j <= n - 1
            A(j+1:m, j + 1) = v(2:m-j+1);
        end
        A
        if j <= n - 2
            [v,b] = house(A(j, j+1:n)');
            % --- added lines ---
            V(j+1:n,j+1:n) = eye(n-j) - b * v * v' * V(j+1:n,j+1:n);
            %---------------------
            A(j:m, j+1:n) = A(j:m, j+1:n) * ( eye(n-j) - b * v * v');
            A(j+1, j+2:n) = v(2:(n-j))';
            disp('---now for the Vs---');
            V
            Aorig * V
            U * Aorig * V
            break;
            disp('here');
        end
        % --- added lines ---
        disp('After 1 iteration');
        B = A
        U
        U * A
        break;
    end
    
    % Case 2: "FAT" MATRIX
else
    last_entry = m;
    for j = 1 : m
        [v,b] = house(A(j:m,j));
        beta_u(j) = b;        
        A(j:m, j:n) = ( eye(m-j+1) - b * (v * v') ) * A(j:m, j:n);
        if j + 1 <= m
            A(j+1:m, j + 1) = v(2:m-j+1);
        end
        if j <= m - 1
            [v,b] = house(A(j, j+1:n)');
            beta_v(j) = b;
            A(j:m, j+1:n) = A(j:m, j+1:n) * ( eye(n-j) - b * v * v');
            A(j+1, j+2:n) = v(2:(n-j))';
        end
    end
    [v,b] = house(A(m, m+1:n)');
    A(m, m+1:n) = A(m, m+1:n) * ( eye(n-m) - b * v * v');
end



% % % Computation of U using backward accumulation
% % k = min(m,n);
% % U = eye(m,m);
% % for j = k : -1 : 1
% %     v = [1; Aorig((j+1):m,j)];
% %     betav = 2/(1 + norm(Aorig((j+1):m,j), 2)^2); % We get the beta's from the stored Householder vectors!
% %     U(j:m,j:m) = U(j:m,j:m) - (betav * (v*v') * U(j:m,j:m));
% % end
% % 
% % U = eye(m,m);
% % for j = k : -1 : 1
% %     v = [1; Aorig((j+1):m,j)];
% %     betav = 2/(1 + norm(Aorig((j+1):m,j), 2)^2); % We get the beta's from the stored Householder vectors!
% %     U(j:m,j:m) = U(j:m,j:m) - (betav * (v*v') * U(j:m,j:m));
% % end
% % 
% % V = 0;



% V = eye(n, n);
% last_entry, m, n
% beta_u
% for j = last_entry :-1 : 1
%     vec = A(j:m, j);
%     disp('loop');
%     U(j:m,j:n) = U(j:m, j:n) - (beta_u(j) * vec ) * (vec' * U(j:m, j:n) ) 
% end
% 
% U
% for j =  n - 2 : -1 : 1
%     vec = A(j, j+1:n)';
%     V(j+1:n,j:n) = V(j:m, j:n) - (beta_v(j) * vec ) * (vec' * V(j+1:n, j:n) ) ;
% end


end