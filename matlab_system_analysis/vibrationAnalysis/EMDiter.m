function [IMF,residue] = EMDiter(x,y)

[U,L]=findPeaks(y);

if length(U)<=1 || length(L)<=1
    IMF = [];
    residue = [];
else

% lb = spline(x(L),y(L),x);
% ub = spline(x(U),y(U),x);
x = reshape(x,1,length(x));
y = reshape(y,1,length(y));

lb = spline([x(1),x(L),x(end)],[y(1),y(L),y(end)],x);
ub = spline([x(1),x(U),x(end)],[y(1),y(U),y(end)],x);

meanLine = (lb+ub)/2;

IMF = y - meanLine;
residue = meanLine;

end
end