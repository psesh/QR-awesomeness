% MODIFIED GRAM-SCHMIDT QR COLUMN PIVOTING ALGORITHM:
%
% Useful references:
% 1. Gene Golub & Charles Van Loan, "Matrix Computations" (2003)
% 2. Aciya Dax, "A modified Gram-schmidt algorithm with iterative
% orthogonalization and column pivoting" (2000)
%
% INPUTS:
%        A: Matrix 
%        fixed_pivots: the first pivots that should not be changed!
%
% OUTPUTS:
%        Q: Orthogonal matrix
%        R: Upper triangular matrix
%        pivots: Column pivots
%
% Copyright (c) 2016 by Pranay Seshadri
%
function [Q,R, pivots] = qr_MGS_pivoting_custom(A, fixed_pivots)

[m,n] = size(A); % Size of "A" -- can set as input!
column_norms = zeros(n,1); % Initialize column norms vector
pivots = 1 : n; % 
r = length(fixed_pivots);

% Create an array with the remainining integers not in fixed_pivots!
combined = [fixed_pivots, pivots];
[~, index] = unique(combined);
pivots = combined(sort(index));
A = A * permutation(pivots');

%---------------------------------------------------------------------
% Step 0:
%---------------------------------------------------------------------
% 1. Compute all the column norms -- this computation is expensive and
% ideal only for the first iteration. [Change me later -- to Pythogras!]
for j = 1 : n
    column_norms(j,1) = norm(A(1:m, j),2)^2;
end

% Now loop!
u = min(m,n);
for k = 1 : u
    
    %---------------------------------------------------------------------
    % Step 0:
    %---------------------------------------------------------------------
    % 2. Find the "j*" column index with the highest column norm
    [~,j_star] = max(column_norms(k:n,1));
    j_star = j_star + (k - 1);

    
    % 3. If j* = k, skip to step 1, else swap columns
    if(k ~= j_star && k > r)
        
        % Swamp columns in A
        temp = A(1:m,k);
        A(1:m,k) = A(1:m,j_star);
        A(1:m,j_star) = temp;
        
        % Swap columns in R
        for i = 1 : k - 1
            temp = R(i,k);
            R(i,k) = R(i,j_star);
            R(i,j_star) = temp;
        end

        % Swap pivots
        temp = pivots(k);
        pivots(k) = pivots(j_star);
        pivots(j_star) = temp;
        
    end
    
    %---------------------------------------------------------------------
    % Step 1: Reorthogonalization
    %---------------------------------------------------------------------
    if( k~=1 )
        for i = 1 : k - 1
            alpha = Q(1:m,i)' * A(1:m,k);
            R(i,k) = R(i,k) + alpha;
            A(1:m,k) = A(1:m,k) - alpha*Q(1:m,i);
        end
    end
    
    %---------------------------------------------------------------------
    % Step 2: Normalization
    %---------------------------------------------------------------------
    R(k,k) =  norm(A(1:m,k),2);
    Q(1:m,k) = A(1:m,k)/R(k,k);
    
    %---------------------------------------------------------------------
    % Step 3: Orthogonalization
    %---------------------------------------------------------------------
    if(k ~= n)
        for j = k + 1 : n
            R(k,j) = Q(1:m,k)' * A(1:m,j);
            A(1:m,j) = A(1:m,j) - R(k,j)* Q(1:m,k);          
            % ---- Seems to be pretty similar to MATLAB -------    
            column_norms(j,1) = norm(A(1:m, j),2)^2;
        end
    end
    
    
end

end
