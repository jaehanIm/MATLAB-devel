function coeff = quadResponseSurface(x,y,data)
% x,y should 'not' be a mesh format

n = length(data(:));
A = zeros(n,6);
B = zeros(n,1);
% c = 1;
% for i = 1:length(x)
%     for j= 1:length(y)
%         A(c,:) = [x(i)^2, y(j)^2, x(i)*y(j), x(i), y(j), 1];
%         B(c) = data(i,j);
%         c = c+1;
%     end
% end

for i = 1:length(x)
    A(i,:) = [x(i)^2, y(i)^2, x(i)*y(i), x(i), y(i), 1];
    B(i) = data(i);
end

coeff = pinv(A) * B;

end