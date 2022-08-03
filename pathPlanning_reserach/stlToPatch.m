function stlToPatch(addr)

tmp = stlread(addr);
FV.faces = tmp.ConnectivityList;
FV.vertices = tmp.Points;

figure(1)
clf
FV = reducepatch(FV,0.5);
p = patch(FV,'FaceColor',       [0.8 0.8 1.0], ...
'EdgeColor',       'none',        ...
'FaceLighting',    'gouraud',     ...
'AmbientStrength', 0.15);
camlight('headlight');
material('shiny')
axis('image');


end