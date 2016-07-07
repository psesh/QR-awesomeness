% Code that computes a step of block QR with column pivoting
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [column_norms, perm] = qr_BlockBlas3QR(m, n, rowk, nb,  A, perm)

% Setup
F(1:n,1:nb) = 0;

% Array for storing column norms
column_norms(1:n) = 0;

% Reduction steps
for j = 1 : nb
    k = rowk + j - 1; % current row index
    
    
    % Compute the column norms
    for r = 1 : n - j + 1
        column_norms(r) = norm(A(k:m, j-1+r), 2);
    end

    % Figure out which index holds the maximum value
    [~,p] = max(column_norms);
    pmax = p + j - 1 % adjust for the location of p
    
    % If the highest norm is 0 then we are done!
    if(column_norms(p) == 0);
        break;
    end
    
    % Below we swap columns of A, elements of the permutation vector
    % and rows of F
    if(p > 1)
        temp = perm(j);
        perm(j) = perm(p);
        perm(p) = temp;
        clear temp;
        
        temp = A(:,j);
        A(:,j) = A(:,pmax);
        A(:,pmax) = temp;
        clear temp;
        
        if(j > 1)
            temp = F(j, 1: j - 1);
            F(j, 1:j-1) = F(pmax, 1 :j-1);
            F(pmax, 1:j-1) = temp;
        end
    end
    
    
    % Update pivot column
    A(k:m,j) = A(k:m,j) - A(k:m, 1:j-1) * F(j, 1:j-1)';
    
    
    % Compute the householder transform -- need to code up this!
    [beta_v, v, tau(j)] = householder_modified(A(k, j), A(k+1:m, j) );
    A(k,j) = 1.0
    if(m - k > 0)
        A(k+1:m, j) = v(2:end);
    end
    
    % Incremental computation of F:
    if(j < n)
        F(j+1:n,j) = tau(j) * A(k:m,j+1:n)' * A(k:m, j);
    end
    
    if(j > 1)
        F(1:n, j) = F(1:n,j) - tau(j) * F(1:n, 1:j-1) * A(k:m, 1:j-1)' * A(k:m,j);
    end
    
    % Updating of pivot row
    A(k,j+1:n) = A(k,j+1:n) - A(k,1:j) * F(j+1:n, 1:j)';
    A(k,j) = beta_v;
    
end
% Block update
if (k < m)
    A(k+1:m, nb+1: n)  =  A(k+1:m, nb+1: n) - A(k+1:m,1:nb)*F(nb+1:n,1:nb)';
end
end

