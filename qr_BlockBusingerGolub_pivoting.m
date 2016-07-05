% Block QR column pivoting using Businger and Golub
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1996)
function [Ablock, column_norms, perm] = qr_BlockBusingerGolub_pivoting(m, n, nb, rowk, A, column_norms, perm)
F(1:n, 1:nb) = 0;
k = 0; % will update in loop!


% Reduction steps
for j = 1 : nb
   k = rowk + j - 1; % current row index
   
   [~,p] = max(column_norms(j:n));
    
    % If the highest norm is 0 then we are done!
    if(column_norms(p) == 0);
        break;
    end
    
    % Swap jth and pth columns!
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
    
    % Update pivot column
    A(k:m,j) = A(k:m,j) - A(k:m,1:j-1) * F(1:j-1,j);
    
    % Reduction
    [v, beta_v] = house(A(k:m, j));
    tau(j) = beta_v;
    Y(j,:) = v;
    Hj = eye(m-k+1) - tau(j) * Y(:,j) * Y(:,j)';
    
    % Incremental computation of F:
    Y
    F(j+1:n,j) = tau(j) * A(j:m,j+1:n)' * Y(j:m,j);
    F(1:n, j) = F(1:n,j) - tau(j) * F(1:n, 1:j-1) * Y(j:m, 1:j-1)' * Y(j:m,j);
    
    % Update of pivot row
    A(k,j+1:n) = A(k,j+1:n) - A(k,1:j) * F(j+1:n, 1:j)';
    
    % Norm downdate
    column_norms(j+1:n) = column_norms(j+1:n) - A(k, j+1: n).^2;
    
end
A(k+1:m, nb+1: n)  = A(k+1:m, nb+1: n) - A(k+1:m,1:nb)*F(nb+1:n,1:nb)';

Ablock = A(k+1:m, nb+1: n) ;
end