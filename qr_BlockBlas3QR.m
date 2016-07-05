% Block QR column pivoting 
%
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1998) SIAM paper
% Coded July 5th 2016
function [Ablock, column_norms, perm] = qr_BlockBlas3QR(m, n, nb, rowk, Asub, perm, column_norms)

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
        
        temp = Asub(:,j);
        Asub(:,j) = Asub(:,p);
        Asub(:,p) = temp;
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
    
    % Update pivot column
    Asub(k:m,j) = Asub(k:m,j) - Asub(k:m, 1:j-1) * F(k, 1:j-1)';
    
    % Reduction
    [v, beta_v] = house(Asub(k:m, j));
    tau(j) = beta_v;
    Y(:,j) = v;
    Hj = eye(m-k+1) - tau(j) * Y(:,j) * Y(:,j)';
    
    % Incremental computation of F:
    F(j+1:n,j) = tau(j) * Asub(j:m,j+1:n)' * Y(j:m,j);
    F(1:n, j) = F(1:n,j) - tau(j) * F(1:n, 1:j-1) * Y(j:m, 1:j-1)' * Y(j:m,j);
    
    % Update of pivot row
    Asub(k,j+1:n) = Asub(k,j+1:n) - Y(k,1:j) * F(j+1:n, 1:j)';
    
    % Norm downdate
    column_norms(j+1:n) = column_norms(j+1:n) - Asub(k, j+1: n).^2;
    
end
Asub(k+1:m, nb+1: n)  = Asub(k+1:m, nb+1: n) - Y(k+1:m,1:nb)*F(nb+1:n,1:nb)';
Ablock = Asub(k+1:m, nb+1: n);
end

