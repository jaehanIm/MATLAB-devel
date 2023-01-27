function coeff = linResponseSurface(x,y,data)
% x,y should 'not' be a mesh format

n = length(data(:));
A = zeros(n,3);
B = zeros(n,1);

for i = 1:length(x)
    A(i,:) = [x(i), y(i), 1];
    B(i) = data(i);
end

coeff = pinv(A) * B;

end
