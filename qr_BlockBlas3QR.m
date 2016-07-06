% Block QR column pivoting 
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [A, column_norms, perm] = qr_BlockBlas3QR(m, n, nb, rowk, A, perm, column_norms)

% Setup
F(1:n,1:nb) = 0;

% Reduction steps
for j = 1 : nb
    k = rowk + j - 1; % current row index
    [~,p] = max(column_norms(j:n));
    
    % If the highest norm is 0 then we are done!
    if(column_norms(p) == 0);
        break;
    end
    
    % Pivoting
    if(j ~= p)
        temp = perm(j);
        perm(j) = perm(p);
        perm(p) = temp;
        clear temp;
        
        temp = A(:,j);
        A(:,j) = A(:,p);
        A(:,p) = temp;
        clear temp;
        
        temp = column_norms(j);
        column_norms(j) = column_norms(p);
        column_norms(p) = temp;
        clear temp;
        
        temp = F(j,:);
        F(j,:) = F(p,:);
        F(p,:) = temp;
        clear temp;
    end
    
    % Reduction
    [v, beta_v] = house(A(k:m, j));
    tau(j) = beta_v;
    Y(:,j) = v;
    
    % Update pivot column
    indices = [k, j, m]
    A(k:m,j) = A(k:m,j-1) - Y(k:m, 1:j-1) * F(j, 1:j-1)';
    disp('got here');

    % Incremental computation of F:
    F(j+1:n,j) = tau(j) * A(j:m,j+1:n)' * Y(j:m,j);
    F(1:n, j) = F(1:n,j) - tau(j) * F(1:n, 1:j-1) * Y(j:m, 1:j-1)' * Y(j:m,j);
    
    % Update of pivot row
    A(k,j+1:n) = A(k,j+1:n) - Y(k,1:j) * F(j+1:n, 1:j)';
    
    % Norm downdate
    column_norms(j+1:n) = column_norms(j+1:n) - A(k, j+1: n).^2;
    
end
A(k+1:m, nb+1: n)  = A(k+1:m, nb+1: n) - Y(k+1:m,1:nb)*F(nb+1:n,1:nb)';

end

