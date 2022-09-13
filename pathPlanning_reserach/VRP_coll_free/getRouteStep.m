function i = getRouteStep(tick, t)

for i = 1:length(tick)
    if tick(i) == 0 && i ~= 1
        i = -1;
        break;
    end
    if t < tick(i)
        i = i-1;
        break;
    end
    if i == length(tick)
        i = -1;
    end
end

end