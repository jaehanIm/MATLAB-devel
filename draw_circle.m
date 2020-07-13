function [x,y,h] = draw_circle(varargin)
% function plotting a circle in your figure
% ** the code does not contain 'hold on' script **
% position is 2x1 vector, r is radius, var contains plot information

pos = varargin{1};
r = varargin{2};
var = varargin{3};
angle = 0:pi/180:2*pi-pi/180;
if length(varargin) == 4
    angle_range = varargin{4};
    angle_range = mod(angle_range,2*pi);
    if angle_range(1) > angle_range(2)
        angle = angle_range(1):-pi/180:angle_range(2);
    else
        angle = angle_range(1):pi/180:angle_range(2);
    end
end

x = pos(1) + r.*cos(angle);
y = pos(2) + r.*sin(angle);

h = plot(x,y,var);

end