function wpList = loadStl(addr,ispDist, reduceFactor)

% tmp = stlread('/home/jaehan/Desktop/F18.stl');
tmp = stlread(addr);

FV.faces = tmp.ConnectivityList;
FV.vertices = tmp.Points;
FV = reducepatch(FV,reduceFactor);

% n = patchnormals(FV);
n = STLVertexNormals(FV.faces,FV.vertices);

figure(99)
clf
patch(FV,'FaceColor',       [0.8 0.8 1.0], ...
'EdgeColor',       'none',        ...
'FaceLighting',    'gouraud',     ...
'AmbientStrength', 0.15);
camlight('headlight');
material('shiny')
axis('image');
alpha(0.5)
hold on
wpList = zeros(size(n,1),3);
for i = 1:size(n,1)
    p1 = FV.vertices(i,:); p2 = FV.vertices(i,:)+ispDist*n(i,:);
%     plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'g-');
    wpList(i,:) = p2;
end

xlabel('x')
ylabel('y')
plot3(wpList(:,1),wpList(:,2),wpList(:,3),'*')

end