stlAddr = '/home/jaehan/Desktop/MATLAB-devel/pathPlanning_reserach/Model/Factory2.stl';

ext1 = [-30:3:25;ones(1,19)*(-17);ones(1,19)*-10]';
ext2 = [ones(1,10)*(-33);-14:3:13;ones(1,10)*(-10)]';
ext3 = [-30:3:0;ones(1,11)*17;ones(1,11)*(-10)]';
ext4 = [3:3:20;ones(1,6)*23;ones(1,6)*(-10)]';
ext5 = [ones(1,11)*26;-13:3:19;ones(1,11)*(-10)]';

extLayer = [ext1;ext2;ext3;ext4;ext5];
extLayer2 = extLayer; extLayer2(:,3) = -6;

ext6 = [-30:3:25;ones(1,19)*17;ones(1,19)*(-2)]';
extLayer3 = [ext1;ext2;ext6;ext5]; extLayer3(:,3) = -2;

ext7 = [ones(1,5)*(-33);0:3:13;ones(1,5)*2]';
extLayer4 = [ext1;ext7;ext6;ext5]; extLayer4(:,3) = 2;
extLayer5 = extLayer4; extLayer5(:,3) = 6;

impNode = loadStl(stlAddr,4);
extLayer6_temp = impNode(impNode(:,3)>10&impNode(:,1)>-28.9&impNode(:,1)<-3&impNode(:,2)<13&impNode(:,2)>-3,:,:);
extLayer6 = extLayer6_temp(((extLayer6_temp(:,1)+20.59).^2+(extLayer6_temp(:,2)-5).^2>49)&((extLayer6_temp(:,1)+11.5).^2+(extLayer6_temp(:,2)-5.3).^2>49),:,:);
selector = zeros(size(extLayer6,1),1);
for i = 1:size(extLayer6,1)
    if rand(1) < 1/12
        selector(i) = 1;
    end
end
extLayer6 = extLayer6(find(selector),:,:); % chimney

extLayer7_temp = impNode(impNode(:,1)<-33 & impNode(:,3)>-5 & impNode(:,2)<0.6 & impNode(:,2)>-12.5,:,:);
extLayer7 = extLayer7_temp(((extLayer7_temp(:,2)+5.85).^2+(extLayer7_temp(:,3)-3.5).^2>25)&((extLayer7_temp(:,2)+5.85).^2+(extLayer7_temp(:,3)-3.5).^2<81),:,:);
selector = zeros(size(extLayer7,1),1);
for i = 1:size(extLayer7,1)
    if rand(1) < 1/15
        selector(i) = 1;
    end
end
extLayer7 = extLayer7(find(selector),:,:); % horizontal chimney

extIO1 = [24,-8,-10;24,-4,-10;24,-8,-6;24,-4,-6];
extIO2 = [24,4,-10;24,7,-10;24,4,-6;24,7,-6];
extIO = [extIO1;extIO2];

extNodes = [extLayer;extLayer2;extLayer3;extLayer4;extLayer5;extLayer6;extLayer7;extIO];

% Interior

intIO1 = [20,-8,-10;20,-4,-10;20,-8,-6;20,-4,-6];
intIO2 = [20,4,-10;20,7,-10;20,4,-6;20,7,-6];
intIO = [intIO1;intIO2];

int1 = [-25:2:19;ones(1,23)*7;ones(1,23)*(-7.5)]'; % corridor
int2 = [3.8+rand(1,25)*16;13.5+rand(1,25)*4.5;ones(1,25)*(-7.5)]'; % small room for 1st floor
int3 = [-27+rand(1,40)*8;-12+rand(1,40)*15;ones(1,40)*(-7.5)]'; % left room for 1st floor
int4 = [-14+rand(1,40)*17.6;-8+rand(1,40)*8;ones(1,40)*(-7.5)]'; % mid room for 1st floor
int5 = [3.6+rand(1,40)*15;-13+rand(1,40)*13;ones(1,40)*(-7.5)]'; % right room for 1st floor

intLayer1 = [int1;int2;int3;int4;int5];
intLayer2 = [int1;int3;int4;int5]; intLayer2(:,3) = 2.5;
stair = [-25*ones(1,6);7*ones(1,6);-6:1.5:1.5]';

intNodes = [intIO;intLayer1;intLayer2;stair];

node = [extNodes;intNodes];

drawStl(stlAddr,1)
hold on
% plot3(node(:,1),node(:,2),node(:,3),'r*')
plot3(extNodes(:,1),extNodes(:,2),extNodes(:,3),'b*')
plot3(intNodes(:,1),intNodes(:,2),intNodes(:,3),'r*')

save('stlnode.mat','node');
