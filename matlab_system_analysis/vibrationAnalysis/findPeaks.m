function [upper, lower] = findPeaks(A)

L = length(A);
prevsgn = 0;
cursgn = 0;
upper = [];
lower = [];

for i = 2:L
    if A(i) > A(i-1)
        cursgn = 1;
    else
        cursgn  = -1;
    end
    
    if prevsgn ~= cursgn
        if prevsgn == 1
            upper = vertcat(upper,i-1);
        elseif prevsgn == -1
            lower = vertcat(lower,i-1);
        end
    end
    
    prevsgn = cursgn;
    
end

end