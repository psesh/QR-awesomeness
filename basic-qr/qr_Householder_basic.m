function [v_store, betav_store] = qr_Householder_basic(A)

% Size of A
[m,n] = size(A);

for j = 1 : n - 1
    [v,betav] = house(A(j:m,j));
    H = (eye(m-j+1) - betav * (v * v') ); % I'd prefer not to use H{j}!
    A(j:m,j:n) =  H * A(j:m,j:n);
    if j < m
        A(j+1:m,j) = v(2:m - j + 1);
    end
    
    % Store
    v_store(:,j) = [zeros(j-1,1) ; v];
    betav_store(j) = betav;
end

end