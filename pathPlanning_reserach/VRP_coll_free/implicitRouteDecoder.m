function route = implicitRouteDecoder(hash)

hash = char(hash);
L = length(hash);
parsePoint = find(hash=='/');
num = length(parsePoint)-1;
route = zeros(1,num);

for i = 1:num
    route(i) = str2double(hash(parsePoint(i)+1:parsePoint(i+1)-1));
end

end