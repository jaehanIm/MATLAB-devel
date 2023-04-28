function hash = implicitRouteEncoder(route)

L = length(route);
hash = char();
for i = 1:L
    hash = horzcat(hash,['/',num2str(route(i))]);  
end
hash = string([hash,'/']);

end