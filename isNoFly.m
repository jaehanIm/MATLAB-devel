function flag = isNoFly(dotpos,wtpos,angle,threshold)
x = dotpos(1); y = dotpos(2);
a = -1/tan(angle);
b = wtpos(2)-a*wtpos(1);

interX = (a*(y-b)+x)/(a^2+1);
interY = interX * a + b;

distance = norm(dotpos - [interX interY]);

if distance <= threshold
    flag = true;
else
    flag = false;
end

end