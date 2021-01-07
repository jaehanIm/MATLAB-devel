%%
% Mission specification calculator for M300
% Lee, Byung-Yoon. NearthLab
% Start date: 200803
% Last modified date: 201104
%

r2d = 180/pi;
d2r = pi/180;

% Camera spec (H20)
% 1/1.7" CMOS, effective 20MP
% Sensor size 7.6 x 5.7 mm, diagonal length 9.5 mm (need to be checked)
% Focal length 6.83 - 119.94 mm (equivalent: 31.7 - 556.2 mm)
% Pixels 5188 x 3888

%%
% FOV and Diagonal focal length
pixelWidth = 5188;
pixelHeight = 3888;
sensorWidth = 7.6; % mm
sensorHeight = 5.7; % mm
sensorRatio = 7.6/5.7;
sensorDiagonalLength = sqrt(sensorWidth^2 + sensorHeight^2);
% focalLength = 6.83;  % mm, min
% focalLength = 119.94;  % mm, max
% focalLength = 103.0;  % mm, target
fovAngle = 2*atan(sensorDiagonalLength/2/focalLength);

% pText = sprintf("For given camera spec (H20)"); disp(pText);
% pText = sprintf(" - Sensor width/height: %f, %f (mm)", sensorWidth, sensorHeight); disp(pText);
% pText = sprintf(" - Pixel width/height: %f, %f", pixelWidth, pixelHeight); disp(pText);
% pText = sprintf(" - Adjustable focal length within 6.83 - 119.94 mm"); disp(pText);
% pText = sprintf(" - Determined focal length: %f (mm)", focalLength); disp(pText);
% pText = sprintf(" - Diagonal FOV angle : %f (deg)\n", (fovAngle*r2d)); disp(pText);

%%
% Inspection distance x m, image spec.
% inspectionDistance = 25.0;
imageDiagonalLength = inspectionDistance*tan(fovAngle);
imageHeight = sqrt(imageDiagonalLength^2/(sensorRatio^2 + 1));
imageWidth = imageHeight * sensorRatio;
imageHeightResolution = imageHeight / pixelHeight * 1000; % mm/pixel
imageWidthResolution = imageWidth / pixelWidth * 1000; % mm/pixel
gsdHeight = (inspectionDistance*1000)*sensorHeight/focalLength/pixelHeight; % mm/pixel
gsdWidth = (inspectionDistance*1000)*sensorWidth/focalLength/pixelWidth; % mm/pixel
imageArea = imageHeight*imageWidth;

% pText = sprintf("For given inspection distance: %f (m)", inspectionDistance); disp(pText);
% pText = sprintf(" - Obtained image width/height: %f / %f (m)", imageWidth, imageHeight); disp(pText);
% pText = sprintf(" - Obtained image width/height resolution: %f / %f (mm/pixel)", ... 
%     imageWidthResolution, imageHeightResolution); disp(pText);
% pText = sprintf(" - Obtained GSD width/height: %f / %f (mm/pixel)", gsdWidth, gsdHeight); disp(pText);
% pText = sprintf(" - Obtained image area: %f (m^2)\n", imageArea); disp(pText);

%%
% Blade spec. (trapezoid)
wtBladeLength = 60;
wtBladeLengthTotal = wtBladeLength*12;
wtLeBladeWidthMin = 1;
wtLeBladeWidthMax = 1;
wtPsBladeWidthMin = 2;
wtPsBladeWidthMax = 5;
wtLeBladeArea = (wtLeBladeWidthMin+wtLeBladeWidthMax)/2*wtBladeLength;
wtPsBladeArea = (wtPsBladeWidthMin+wtPsBladeWidthMax)/2*wtBladeLength;
wtBladeArea = wtLeBladeArea*2 + wtPsBladeArea*2;
wtArea_6side = wtBladeArea*3/2;
wtArea_12side = wtBladeArea*3;

% pText = sprintf("For given WT blade length: %f (m)", wtBladeLength); disp(pText);
% pText = sprintf("For given WT blade spec LE width min/max: %f / %f (m)", ...
%     wtLeBladeWidthMin, wtLeBladeWidthMax); disp(pText);
% pText = sprintf("For given WT blade spec PS width min/max: %f / %f (m)", ...
%     wtPsBladeWidthMin, wtPsBladeWidthMax); disp(pText);
% pText = sprintf(" - Total blade length: %f (m)", wtBladeLengthTotal); disp(pText);
% pText = sprintf(" - Area of one blade: %f (m^2)", wtBladeArea); disp(pText);
% pText = sprintf(" - Area of 12-side blade: %f (m^2)", wtArea_12side); disp(pText);
% pText = sprintf(" - Area of 6-side blade: %f (m^2)\n", wtArea_12side); disp(pText);

%%
% Estimated shoot time
reqShootNumber = wtArea_6side/imageArea;
shootNumberMargin = 2.0;
reqShootNumber = reqShootNumber * shootNumberMargin;
shootTime = 1.5; % sec
reqShootTime = reqShootNumber * shootTime;
reqMissionTimeMinute = reqShootTime / 60;

% pText = sprintf("For given shoot time: %f (sec)", shootTime); disp(pText);
% pText = sprintf(" - Required shoot number: %f", reqShootNumber); disp(pText);
% pText = sprintf(" - Required inspection time: %f (min)\n", reqMissionTimeMinute); disp(pText);

%%
% Gimbal control range
% Pan +-320 deg, tilt -120 ~ +30 deg
gimbalTiltAngle = 13.0 * d2r;
imageTiltAngle = gimbalTiltAngle + fovAngle/2;
distFromGimbalZero = inspectionDistance*tan(gimbalTiltAngle);
imageCoverageFromGimbalZero = distFromGimbalZero + imageHeight/2;
distanceBwShootPoint = imageCoverageFromGimbalZero*2;
shootPointNumber = wtBladeLengthTotal/distanceBwShootPoint;
shootPointNumberPerBladeOneSide = shootPointNumber/12;
shootNumberPerShootPoint = reqShootNumber/shootPointNumber;

% pText = sprintf("For acceptable gimbal/image tilt angle against blade: %f / %f (deg)", ...
%     gimbalTiltAngle*r2d, imageTiltAngle*r2d); disp(pText);
% pText = sprintf(" - Image coverage distance: +- %f (m)", imageCoverageFromGimbalZero); disp(pText);
% pText = sprintf(" - Min. distance b/w shoot points: %f (m)", distanceBwShootPoint); disp(pText);
% pText = sprintf(" - Number of shoot points: %f", shootPointNumber); disp(pText);
% pText = sprintf(" - Shoot point # per blade one side: %f", shootPointNumberPerBladeOneSide); disp(pText);
% pText = sprintf(" - Shoot number per each shoot point: %f\n", shootNumberPerShootPoint); disp(pText);

%%
% Flight Time Estimation for 1 sortie mission(6-side)
% pText = sprintf("Flight Time Estimation"); disp(pText);

% WT modeling time(WT nose)
reqWtNoseModelTime = 30.0;
% pText = sprintf(" - Required WT nose modelling time: %f (sec)", reqWtNoseModelTime); disp(pText);

% WT modeling time(WT blade)
reqWtBladeModelTimeOnce = 5.0;
reqWtBladeModelTime = reqWtBladeModelTimeOnce * ceil(shootPointNumberPerBladeOneSide) * 6;
% pText = sprintf(" - Required WT blade modelling time: %f (sec)", reqWtBladeModelTime); disp(pText);

% Shoot time
reqShootTime;
% pText = sprintf(" - Required inspection time: %f (sec)", reqShootTime); disp(pText);

% Hopping time
vehicleHoppingVelocity = 2.0; % m/s
vehicleHoppingDistance = distanceBwShootPoint*(ceil(shootPointNumberPerBladeOneSide)-1)*6;
reqHoppingTime = vehicleHoppingDistance/vehicleHoppingVelocity;
% pText = sprintf(" - Required hopping time: %f (sec)", reqHoppingTime); disp(pText);

% Surface transition time
vehicleSurfaceTransitionVelocity = 4.0; % m/s
surfaceTransitionDistanceOnce = 2*pi*inspectionDistance;
surfaceTransitionDistance = surfaceTransitionDistanceOnce * 6;
reqSurfaceTransitionTime = surfaceTransitionDistance/vehicleSurfaceTransitionVelocity;
% pText = sprintf(" - Required surface transition time: %f (sec)", reqSurfaceTransitionTime); disp(pText);

% Take-off/Landing time
reqToLgTime = 180.0;
% pText = sprintf(" - Required T/O, L/G time: %f (sec)", reqToLgTime); disp(pText);

% Total flight time
reqTotalFlightTime = reqWtNoseModelTime + reqWtBladeModelTime + reqShootTime ...
    + reqHoppingTime + reqSurfaceTransitionTime + reqSurfaceTransitionTime;
reqTotalFlightTimeMin = reqTotalFlightTime/60.0;
% pText = sprintf(" - Required Total Flight Time: %f (sec)", reqTotalFlightTime); disp(pText);
% pText = sprintf(" - Required Total Flight Time: %f (min)", reqTotalFlightTimeMin); disp(pText);




