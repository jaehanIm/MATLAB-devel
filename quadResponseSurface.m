function coeff = quadResponseSurface(x,y,data)

n = length(data(:));
A = zeros(n,6);
B = zeros(n,1);
c = 1;
for i = 1:length(x)
    for j= 1:length(y)
        A(c,:) = [x(i)^2, y(j)^2, x(i)*y(j), x(i), y(j), 1];
        B(c) = data(i,j);
        c = c+1;
    end
end

coeff = pinv(A) * B;

end