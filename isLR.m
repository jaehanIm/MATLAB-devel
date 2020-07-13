function [flag, distance] = isLR(dotpos,wtpos,angle)
% position of 'dotpos' relative to WT yaw angle line
x = dotpos(1); y = dotpos(2);
angle = mod(angle, 2*pi);
a = tan(angle);
b = wtpos(2)-a*wtpos(1);

interX = (a*(y-b)+x)/(a^2+1);
interY = interX * a + b;

out = cross([x y 0]-[interX interY 0],[cos(angle) sin(angle) 0]);
distance = norm(dotpos - [interX interY]);

if out(3) > 0
    flag = "right";
else
    flag = "left";
end

end