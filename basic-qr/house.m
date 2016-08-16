% Computation of the Householder vector. For an input x in R^(m x 1), t
% this function computes v in R^(m x 1) with v(1) = 1 and a constant betav 
% such that P = eye(m) - betav*v*v' is orthogonal. 
%
% Copyright (c) 2016 by Pranay Seshadri
function [v,betav] = house(x)

m = length(x);
sigmav = x(2:m)' * x(2:m);
v = [1; x(2:m) ];

if(length(x) == 1)
    v = 1;
    betav = x(1);
elseif(sigmav == 0 && x(1) >= 0 )
    betav = 0;
elseif(sigmav == 0 && x(1) < 0 )
    betav = -2;
else
    muv = sqrt(x(1)^2 + sigmav);
    if x(1) <= 0
        v(1) = x(1) - muv;
    else
        v(1) = -sigmav/(x(1) + muv) ;
    end
    betav = 2*v(1)^2 / (sigmav + v(1)^2);
    v = v/v(1);
end
