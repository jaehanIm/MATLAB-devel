addr = '/home/jaehan/nearthlab/niv_sim/nlab_gazebo/models/wt_model_2/meshes/Industrial_031.stl';

tmp = stlread(addr);
FV.faces = tmp.ConnectivityList;
FV.vertices = tmp.Points;

scale = 0.1;
rotM = angle2dcm(pi/2,0,0,'XYZ')';
FV.vertices = rotM * FV.vertices';
FV.vertices = FV.vertices';
FV.vertices = FV.vertices*scale;
FV.vertices = FV.vertices + [0,+37,0];

figure(12)
hold on
patch(FV,'FaceColor',       [0.8 0.8 1.0], ...
'EdgeColor',       'none',        ...
'FaceLighting',    'gouraud',     ...
'AmbientStrength', 0.15);
camlight('headlight');
material('dull')
axis('image');
alpha(0.7)
xlabel('x')
ylabel('y')
view(30,40)