function [angle, direc] = anglerel(t1, t2)
t1 = mod(t1,2*pi); t2 = mod(t2,2*pi);

if t2 - t1 >= pi
    direc = "right";
    angle = t2 - t1 - 2 * pi;
elseif t2 - t1 >= -pi
    angle = t2 - t1;
    if t2 - t1 >= 0
        direc = "left";
    else
        direc = "right";
    end
else
    direc = "left";
    angle = t2 - t1 + 2 * pi;
end

if t2 == t1
    angle = 0;
    direc = "identical";
end

end