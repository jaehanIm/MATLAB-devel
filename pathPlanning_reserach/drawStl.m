function drawStl(addr,figNum)

tmp = stlread(addr);
FV.faces = tmp.ConnectivityList;
FV.vertices = tmp.Points;

figure(figNum)
patch(FV,'FaceColor',       [0.8 0.8 1.0], ...
'EdgeColor',       'none',        ...
'FaceLighting',    'gouraud',     ...
'AmbientStrength', 0.15);
camlight('headlight');
material('shiny')
axis('image');
alpha(0.8)

end