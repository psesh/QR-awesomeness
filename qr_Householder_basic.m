function H = qr_Householder_basic(A)

% Size of A
[m,n] = size(A);

for j = 1 : n - 1
    [v,betav] = house(A(j:m,j));
    H{j} = (eye(m-j+1) - betav * (v * v') ); % I'd prefer not to use H{j}!
    A(j:m,j:n) =  H{j} * A(j:m,j:n);
    if j < m
        A(j+1:m,j) = v(2:m - j + 1);
    end
end

end