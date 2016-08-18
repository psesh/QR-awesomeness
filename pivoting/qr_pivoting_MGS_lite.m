% MODIFIED GRAM-SCHMIDT QR COLUMN PIVOTING ALGORITHM:
% See:
% 1. Gene Golub & Charles Van Loan, "Matrix Computations" (2003)
% 2. Aciya Dax, "A modified Gram-schmidt algorithm with iterative
% orthogonalization and column pivoting" (2000)
%
% Copyright (c) 2016 by Pranay Seshadri
%
function pivots = qr_pivoting_MGS_lite(A)

% preliminaries
[m,n] = size(A); 
column_norms = zeros(n,1);
pivots = 1 : n;
u = min(m,n);

% computation of column norms
for j = 1 : n
    column_norms(j,1) = norm(A(1:m, j),2)^2;
end


for k = 1 : u
    
    % compute highest norm
    [~,j_star] = max(column_norms(k:n,1));
    j_star = j_star + (k - 1);

    % Retrieve the k-th column of A
    a_k = A(:,k);
    
    % swap
    if(k ~= j_star)
        
        % collect a_jstar
        a_jstar = A(:,j_star);
        
        % Swap
        temp = a_k;
        a_k = a_jstar;
        a_jstar = temp;

        temp = pivots(k);
        pivots(k) = pivots(j_star);
        pivots(j_star) = temp;
        
    end
    
    A(:,k) = a_k;
    A(:, j_star) = a_jstar;
    
    % orthognalization
    if(k ~= n)
        for j = k + 1 : n
            a_j = A(:,j);
            a_j = a_j -  (a_k/norm(a_k,2))' * a_j * (a_k/norm(a_k,2));  
            column_norms(j,1) = norm(a_j,2)^2;
            A(:,j) = a_j;
        end
    end 
    
    % reorthogonalization
    if( k~=1 )
        for i = 1 : k - 1
            a_i = A(:,i);
            a_k = a_k - (a_i/norm(a_i,2))' * a_k * (a_i/norm(a_i,2));
        end
    end
   
    A(:,k) = a_k;
    
end

end
