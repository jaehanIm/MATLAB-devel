%% Parameter setting

px = 6 * fovFactor; % fov x length
py = 4 * fovFactor;
ov_x = 0.1; % overlap length
ov_y = 0.1;
voxel_width = 0.2; % voxel distance

eff_x = px - ov_x; % Effective fov size
eff_y = py - ov_y;

%% Load data (Inspection Area Loading)
area = imread('area.png');
area = area~=0;
area = ones(size(area));

area_xl = size(area,2)*voxel_width; % Acutal size[m] of area
area_yl = size(area,1)*voxel_width;
[voxelPosX, voxelPosY] = meshgrid(linspace(voxel_width, voxel_width*(size(area,2)+1),size(area,2)),linspace(voxel_width, voxel_width*(size(area,1)+1),size(area,1)));
areaPosX = voxelPosX;
areaPosY = voxelPosY;

%% Grid Initialization
gridnum_x = ceil(area_xl/eff_x); % required grid size
gridnum_y = ceil(area_yl/eff_y);

gridPosX = []; gridPosY = []; gridValue = [];

[gridPosX, gridPosY] = meshgrid(eff_x/2:eff_x:eff_x*gridnum_x-eff_x/2,eff_y/2:eff_y:eff_y*gridnum_y-eff_y/2); % grid coordinate
gridEdgeX = [gridPosX(1,:)-eff_x/2,gridPosX(end)+eff_x/2];
gridEdgeY = [gridPosY(:,1)'-eff_y/2,gridPosY(end)+eff_y/2];

%% Curvature generation (Voxel Data Loading)
voxelData = peaks(length(area))*mapheight;
voxelData = voxelData(1:size(area,1),1:size(area,2));
voxelFilterData = area.*voxelData;

gridValue = griddata(voxelPosX(:),voxelPosY(:),voxelData(:),gridPosX,gridPosY);

voxelFilterData(area==0) = nan; % neglect unselected area

%% Voxel to grid localization
vg = {};
for i = 1:size(gridPosX,2)
    for j = 1:size(gridPosX,1)
        vg(i,j).x = voxelPosX(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelData));
        vg(i,j).y = voxelPosY(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelData));
        vg(i,j).data = voxelData(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelData));
        vg(i,j).filterData = voxelData(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelFilterData));
        vg(i,j).filterX = voxelPosX(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelFilterData));
        vg(i,j).filterY = voxelPosY(voxelPosX > gridEdgeX(i) & voxelPosX < gridEdgeX(i+1) & voxelPosY > gridEdgeY(j) & voxelPosY < gridEdgeY(j+1) & ~isnan(voxelFilterData));
        vg(i,j).null = 0;
        if isempty(vg(i,j).filterData)
           vg(i,j).null = 1;
        end
    end
end

%% Local RSM gen.
gridNormVector = [];
gridNormVector_lin = [];
gridNormVector{size(gridValue,1),size(gridValue,2)} = [];
gridNormVector_lin{size(gridValue,1),size(gridValue,2)} = [];

for i = 1:size(vg,1)
    for j = 1:size(vg,2)
        if vg(i,j).null == 0
            % Linear approximation
            coeff = linResponseSurface(vg(i,j).x, vg(i,j).y, vg(i,j).data);
            gridNormVector_lin{j,i} = [coeff(1),coeff(2),-1];

            % Quadratic approximation
            coeff = quadResponseSurface(vg(i,j).x, vg(i,j).y, vg(i,j).data);
            temp_x = gridPosX(1,i);
            temp_y = gridPosY(j,1);
            gridNormVector{j,i} = slopeQuadResponseSurf(coeff,temp_x,temp_y);
        else
            gridNormVector{j,i} = ones(1,3)*nan;
            gridNormVector_lin{j,i} = ones(1,3)*nan;
        end
    end
end

%% Air point gen.
airPosX = zeros(size(gridValue)); airPosY = zeros(size(gridValue)); airPosZ = zeros(size(gridValue));
airPosX_lin = zeros(size(gridValue)); airPosY_lin = zeros(size(gridValue)); airPosZ_lin = zeros(size(gridValue));
for i = 1:size(gridNormVector,1)
    for j = 1:size(gridNormVector,2)
        airPosX(i,j) = gridPosX(i,j) - inpection_dist * gridNormVector{i,j}(1);
        airPosY(i,j) = gridPosY(i,j) - inpection_dist * gridNormVector{i,j}(2);
        airPosZ(i,j) = gridValue(i,j) - inpection_dist * gridNormVector{i,j}(3);
        airPosX_lin(i,j) = gridPosX(i,j) - inpection_dist * gridNormVector_lin{i,j}(1);
        airPosY_lin(i,j) = gridPosY(i,j) - inpection_dist * gridNormVector_lin{i,j}(2);
        airPosZ_lin(i,j) = gridValue(i,j) - inpection_dist * gridNormVector_lin{i,j}(3);
    end
end
