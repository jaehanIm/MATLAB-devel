function fix = nodeInterpreter(num)

if num == 1 || num == 0
    home = 1;
    rt = 0;
    blNum = 0;
    fixNum = 0;
else
    index = num-1;
    home = 0;
    blNum = ceil(index/8);
    leftover = index - (blNum-1)*8;
    if leftover > 4
        rt = 2;
        fixNum = leftover - 4;
    else
        rt = 1;
        fixNum = leftover;
    end
end

fix.home = home;
fix.rt = rt;
fix.blNum = blNum;
fix.fixNum = fixNum;

end