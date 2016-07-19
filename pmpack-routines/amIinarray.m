% Simple function that takes in a value and an array and determines
% if the value is in the array.
%
% Copyright (c) 2016 by Pranay Seshadri
% University of Cambridge
%
function [yesOrNo, location] = amIinarray(value, array)

for i = 1 : length(array)
    if value == array(i)
        location = i;
        yesOrNo = 'true';
        break;
    else
        location = 0;
        yesOrNo = 'false';
    end
end

end