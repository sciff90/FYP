function err = funcs(y1,y2,num_samples)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
err = norm(1/num_samples.*sum(y1-y2).^2)

end

