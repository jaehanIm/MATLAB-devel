function coeff = linResponseSurface(x,y,data)
% x,y should 'not' be a mesh format

n = length(data(:));
A = zeros(n,3);
B = zeros(n,1);
c = 1;
for i = 1:length(x)
    for j= 1:length(y)
        A(c,:) = [x(i), y(j), 1];
        B(c) = data(i,j);
        c = c+1;
    end
end

coeff = pinv(A) * B;

end