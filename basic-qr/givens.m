% Algorithm to compute "c" and "s" for a Givens Rotation
% Copyright (c) 2016 by Pranay Seshadri
function [c,s] = givens(a,b)
if(b==0)
    c = 1;
    s = 0;
else
    if(abs(b) > abs(a) )
        tau = -a/b;
        s = 1/sqrt(1 + tau^2);
        c = s * tau;
    else
        tau = -b/a;
        c = 1/sqrt(1 + tau^2);
        s = c * tau;
    end
end

end
