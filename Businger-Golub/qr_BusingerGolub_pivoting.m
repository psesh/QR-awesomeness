% Variant of Businger and Golub's QR with column pivoting
% Copyright (c) 2016 by Pranay Seshadri
%
% See Quintana-Orti's et al. (1996)
function [R,perm] = qr_BusingerGolub_pivoting(A)
[m,n] = size(A);
perm = 1 : n;

% Compute column norms
for j = 1 : n
   column_norms(j) = norm(A(:,j),2); 
end

% Reduction steps
for j = 1 : n - 1
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
        
    end
    
    % Reduction -- compute Householder matrix
    [v, beta_v] = house(A(j:m, j));
    Hj = (eye(m-j+1) - beta_v * (v * v') );
    A(j:m,j+1:n) = Hj * A(j:m,j+1:n);
     
    
    % Norm downdate
    column_norms(j+1:n) = column_norms(j+1:n) - A(j,j+1:n).^2;
   
end
R = triu(A);

end