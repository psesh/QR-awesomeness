function output = funceval(f, pts)
% Evaluates a multivariate function give multivariate points:
% 
% Pranay Seshadri
% Copyright 2012
% ps583@cam.ac.uk
%
% Sample input:
% f = @(X) exp(X(1) + X(2));
% x1  = -1:0.1:1; x2  = -1:0.1:1;
% points = [x1;x2]';
% output = funceval(f, points);
%
%
[row,~] = size(pts);
output = zeros(row, 1);
for i = 1 : row % For the number of points in each direction
    output(i) = feval(f, pts(i,:));
end
end