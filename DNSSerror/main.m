clc
% d2r = pi/180;
% disp(['AP _ 103712 _ 1before'])
% dcm1 = angle2dcm(70.00565*d2r,-22.48895*d2r,94.90314*d2r,'ZXY');
% dcm2 = angle2dcm(67.68523*d2r,-24.77648*d2r,85.79105*d2r,'ZXY');
% refVec = [0,0,1]';
% rootToTipVec = [0.166671 0.904287 0.393048]';
% surfaceNormal = [-0.553416 -0.244122 0.796326]';
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% projRoot = rootToTipVec - surfaceNormal * dot(rootToTipVec, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),projRoot/norm(projRoot))
% dot(projVec2/norm(projVec2),projRoot/norm(projRoot))])/d2r
% disp(['AP _ 103712 _ 1after'])
% dcm1 = angle2dcm(-83.67928*d2r,21.16283*d2r,79.74955*d2r,'ZXY');
% dcm2 = angle2dcm(-82.8472*d2r,22.40572*d2r,66.66727*d2r,'ZXY');
% refVec = [0,0,1]';
% rootToTipVec = [-0.145034 -0.786894 0.599802]';
% surfaceNormal = [-0.397564 0.601469 0.692949]';
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% disp(['AP _ 103712 _ 2before'])
% dcm1 = angle2dcm(69.08791*d2r,-23.11333*d2r,94.90314*d2r,'ZXY');
% dcm2 = angle2dcm(67.49952*d2r,-24.38341*d2r,85.79105*d2r,'ZXY');
% dcm3 = angle2dcm(66.31706*d2r,-25.55473*d2r,74.87302*d2r,'ZXY');
% dcm4 = angle2dcm((64.28296-97.5233)*d2r,66.4003*d2r,-101.324*d2r,'ZYX');
% refVec = [0,0,1]';
% rootToTipVec = [0.166291 0.897633 0.40817]';
% surfaceNormal = [-0.556024 -0.256494 0.7906]'; %NEW
% surfaceNormal = [-0.326306 -0.340521 0.881799]'; %OLD
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% imgVec3 = dcm3'*refVec;
% imgVec4 = dcm4'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% projVec3 = imgVec3 - surfaceNormal * dot(imgVec3, surfaceNormal);
% projVec4 = imgVec4 - surfaceNormal * dot(imgVec4, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)
% dot(projVec3/norm(projVec3),rootToTipVec)
% dot(projVec4/norm(projVec4),rootToTipVec)])./d2r
% disp(['AP _ 103712 _ 2after'])
% dcm1 = angle2dcm(-83.56999*d2r,21.82076*d2r,80.4295*d2r,'ZXY');
% dcm2 = angle2dcm(-82.93778*d2r,22.85433*d2r,67.4429*d2r,'ZXY');
% refVec = [0,0,1]';
% rootToTipVec = [-0.147535 -0.796388 0.586515]';
% surfaceNormal = [-0.399111 0.590516 0.701427]';
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% % SILS _ 230103_175114
% disp(['SILS _ 230103_175114'])
% dcm1 = angle2dcm(117.15*d2r+0.0694164,-0.392352,1.53681,'ZXY');
% dcm2 = angle2dcm(117.15*d2r+0.0479883,-0.418067,1.30041,'ZXY');
% rootToTipVec = [-0.614852 0.615178 0.493471]';
% surfaceNormal = [-0.0515393 -0.655728 0.753236]';
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% % AG _ incheon _ 135611
% disp(['AG _ incheon _ 135611'])
% rootToTipVec = [0.860937 -0.0464298 0.506588]';
% surfaceNormal = [-0.465007 -0.475649 0.746677]';
% dcm1 = angle2dcm(12.02899*d2r,23.58947*d2r,73.57244*d2r,'ZXY');
% dcm2 = angle2dcm(14.98703*d2r,24.4325*d2r,55.09671*d2r,'ZXY');
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% % AG 105416 _ 2before
% rootToTipVec = [0.0943585 0.909571 0.404695]';
% surfaceNormal = [-0.533495 -0.297012 0.791938]';
% dcm1 = angle2dcm(72.47069*d2r,-23.9483*d2r,86.07281*d2r,'ZXY');
% dcm2 = angle2dcm(71.17461*d2r,-25.36281*d2r,75.16406*d2r,'ZXY');
% dcm3 = angle2dcm(68.96107*d2r,-26.29431*d2r,62.04756*d2r,'ZXY');
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% imgVec3 = dcm3'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% projVec3 = imgVec3 - surfaceNormal * dot(imgVec3, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)
% dot(projVec3/norm(projVec3),rootToTipVec)])/d2r
% % AG 105416 _ 2after
% rootToTipVec = [-0.0833435 -0.803391 0.589591]';
% surfaceNormal = [-0.444644 0.559468 0.69949]';
% dcm1 = angle2dcm(-79.31422*d2r,21.39286*d2r,80.37103*d2r,'ZXY');
% dcm2 = angle2dcm(-78.51737*d2r,22.61332*d2r,67.31461*d2r,'ZXY');
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% disp(['afterfix test fixed'])
% rootToTipVec = [-0.00085347 0.869666 0.49364]';
% surfaceNormal = [-0.486809 -0.447186 0.750361]';
% dcm1 = angle2dcm(76.4227*d2r,-22.497*d2r,88.0789*d2r,'ZXY');
% dcm2 = angle2dcm(75.18592*d2r,-23.9847*d2r,74.5279*d2r,'ZXY');
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
% disp(['afterfix test unfixed'])
% rootToTipVec = [-0.0259822 0.866053 0.499277]';
% surfaceNormal = [-0.49958 -0.427995 0.753153]';
% dcm1 = angle2dcm(77.80559*d2r,-22.398*d2r,87.6383*d2r,'ZXY');
% dcm2 = angle2dcm(76.61921*d2r,-23.8666*d2r,74.1096*d2r,'ZXY');
% imgVec1 = dcm1'*refVec;
% imgVec2 = dcm2'*refVec;
% projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
% projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
% acos([dot(projVec1/norm(projVec1),rootToTipVec)
% dot(projVec2/norm(projVec2),rootToTipVec)])/d2r
disp(['SILS test'])
dcm1 = angle2dcm(3.95734*d2r,-22.398*d2r,87.6383*d2r,'ZXY'); %old
% dcm2 = angle2dcm(2.27375*d2r,-22.497*d2r,88.0789*d2r,'ZXY'); %new
dcm2 = angle2dcm(2.28715*d2r,-22.4756*d2r,88.0768*d2r,'ZXY'); %new
% dcm2 = angle2dcm(3.979*d2r,-22.48*d2r,88.07*d2r,'ZXY'); %old
refVec = [0,0,1]';
surfaceNormal = [-0.486809 -0.447186 0.750361]'; %NEW
imgVec1 = dcm1'*refVec;
imgVec2 = dcm2'*refVec;
projVec1 = imgVec1 - surfaceNormal * dot(imgVec1, surfaceNormal);
projVec2 = imgVec2 - surfaceNormal * dot(imgVec2, surfaceNormal);
acos(dot(projVec1/norm(projVec1), projVec2/norm(projVec2)))/d2r
normalize(cross(projVec1,projVec2),'norm')

pt1 = dcm1'*[1,0,0]'*9;
pt2 = dcm2'*[1,0,0]'*9;
pt1-pt2

disp(['rotation test'])
v1 = [14.99937, 12.59579, 5.06541]'; %우측 
v2 = [15.38097, -7.434131, 5.080467]'; %좌측
dcm1 = angle2dcm(-6.85573*d2r,0,0,'XYZ');
v1 = dcm1*v1
v2 = dcm1*v2
diffvec = v1-v2;
atan2(diffvec(1),diffvec(2))/d2r

v1 = [14.99937, 12.59579, 5.06541]'; %우측 
temp = v1(1:2)-[16.0033,-0.6827838]';
atan2(temp(1),temp(2))/d2r
%% damn it
% curYawTarget = 0 * d2r;
% vehPos = [200.4409, -140.1205, -23.82028]';
% shootNormalW = [-0.151293, 0.604197, 0.782341]';
% shootNormalN = [-0.399111, 0.590516, 0.701427]';
% shootPos = [204.519 -143.452 -34.5901]';
% posToTarget = shootPos - vehPos;
% reqAngle = calculator(vehPos, shootNormalW, posToTarget,curYawTarget);
% reqAngle/d2r
% 
% reqDcm = angle2dcm(reqAngle(3),reqAngle(2),reqAngle(1),'ZYX');
% [yaw, roll, pitch] = dcm2angle(reqDcm,'ZXY');
% gimbalAngle = [roll,pitch,yaw]/d2r
% 
% %test code
% rootToTipVec = [0.166291 0.897633 0.40817]';
% dcm1 = angle2dcm(yaw+curYawTarget,roll,pitch,'ZXY');
% imgVec1 = dcm1'*[0,0,1]';
% projVec1 = imgVec1 - shootNormalW * dot(imgVec1, shootNormalW);
% dot(projVec1/norm(projVec1),rootToTipVec)
% 
% function out = calculator(vehPos, shootNormal, posToTarget,curYawTarget)
% d2r = pi/180;
% wtHdg = -10.4429 * d2r;
% blRot = (114.09-90) * d2r;
% dcmI2ShootSurface = -angle2dcm(wtHdg,blRot,60*d2r,'ZXY'); %75 wrong %60 new
% shootNormal = dcmI2ShootSurface'*[1,0,0]'
% rootToTipVec = [0.166291 0.897633 0.40817]';
% 
% posToTarget0 = posToTarget(1:2);
% shootYawAngleReq = wrapToPi(atan2(posToTarget(2),posToTarget(1))-curYawTarget);
% shootPitchAngleReq = acos(norm(posToTarget0)/norm(posToTarget));
% 
% dcmI2FovPV = angle2dcm(shootYawAngleReq+curYawTarget,shootPitchAngleReq,0,'ZYX');
% shootYPNormalVecNed = dcmI2FovPV' * [-1.0, 0.0, 0.0]';
% surfaceVecDotProduct = dot(shootYPNormalVecNed, shootNormal);
% surfaceTiltAngle = acos(surfaceVecDotProduct);
% 
% fovRefVecNed = dcmI2FovPV' * [0,0,1]';
% fovRefVecShootSurface = -dcmI2ShootSurface*fovRefVecNed;
% fovRefVecShootSurface(3) = fovRefVecShootSurface(3);
% revVecShootSurfaceFlat = fovRefVecShootSurface;
% revVecShootSurfaceFlat(1) = 0;
% 
% rootToTipVecShootSurface = -dcmI2ShootSurface * rootToTipVec;
% rootToTipVecShootSurface(3) = -rootToTipVecShootSurface(3);
% rootToTipVecShootSurfaceFlat = rootToTipVecShootSurface;
% rootToTipVecShootSurfaceFlat(1) = 0;
% 
% projectedVecAngleDiff = abs(acos(dot(rootToTipVecShootSurfaceFlat,revVecShootSurfaceFlat)/norm(rootToTipVecShootSurfaceFlat)/norm(revVecShootSurfaceFlat)));
% shootRollAngleReq_ = atan(tan(projectedVecAngleDiff)*cos(surfaceTiltAngle));
% if (projectedVecAngleDiff > pi/2)
%     shootRollAngleReq_ = shootRollAngleReq_ + sign(projectedVecAngleDiff) * pi;
% end
% 
% rotationDirection = -sign(cross(rootToTipVecShootSurfaceFlat,revVecShootSurfaceFlat));
% shootRollAngleReq_ = rotationDirection * shootRollAngleReq_;
% shootRollAngleReq_(2:end) = [];
% 
% out = [shootRollAngleReq_,shootPitchAngleReq,shootYawAngleReq];
% end