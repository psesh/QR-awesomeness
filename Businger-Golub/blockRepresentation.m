function [W, Y] = blockRepresentation(v, betav)
% Algorithm computes matrices W, Y in R^{m, n} such that Q = Im - WY'
Y = v(:,1);
W = betav(1) * v(:,1);
[rows, cols] = size(v);

for j = 2 : cols
   M = W * Y';
   [r2,c2] = size(M);
   z = betav(j) * (eye(r2,c2) - M) * v(:,j);
   W = [W, z];
   Y = [Y, v(:,j)];
end


end