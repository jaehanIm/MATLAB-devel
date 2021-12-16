function out = inBetween(series,reference)

L = length(series);
prevsgn = sign(series(1)-reference);
for i = 2:L
    cursgn = sign(series(i)-reference);
    if cursgn * prevsgn == -1
        out = [i-1, i];
        break;
    elseif cursgn*prevsgn == 0
        out = [i i];
        break;
    end
    prevsgn = cursgn;
end


end