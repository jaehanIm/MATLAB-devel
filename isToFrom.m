function [flag, distance] = isToFrom(dotpos,wtpos,angle)
x = dotpos(1); y = dotpos(2);
a = -1/tan(angle);
b = wtpos(2)-a*wtpos(1);

interX = (a*(y-b)+x)/(a^2+1);
interY = interX * a + b;

out = dot([x y]-[interX interY],[cos(angle) sin(angle)]);
distance = norm(dotpos - [interX interY]);

if out > 0
    % coincide w/ WT angle
    flag = "front";
else
    flag = "back";
end

end