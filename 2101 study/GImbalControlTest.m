%% Parameter setting
d2r = pi/180; r2d = 1/d2r;

yp = -25.7723; yp = yp * d2r; % Plane pitch angle
zp = -31.6616; zp = zp * d2r; % Plane yaw angle
yr = -32.0866; yr = yr*d2r; % Required pitch angle
zr = -18.2713; zr = zr*d2r; % Required yaw angle

n = 100; % Sample size

%% Calculate gimbal attitude

if yr == 0 % Special case that there is no pitch request -> dead zone design
    xg = zeros(1,n);
    yg = zeros(1,n);
    zg = zr * ones(1,n);
else
    % Generate x, y, z samples
    zg = linspace(-88,88,n) * d2r;
%     zg = linspace(zr*r2d-10,zr*r2d+10,n) * d2r;
%     zg = -17.5642 * pi/180;
    dz = zr - zg;
    alpha = cos(yr).*cos(dz);
    beta = cos(yr).*sin(dz);
    xg = asin(sqrt(beta.^2./(1-alpha.^2))) .*dz ./abs(dz) .*yr ./abs(yr);
    yg = acos(alpha) * yr / abs(yr);
    
    % Calculate roll directional condition
    cxg = cos(xg);
    syp = sin(yg);
    invOAF = -sin(zg-zp);
    sxg = sin(xg);
    cyp = cos(yg);
    
    direcCond = cxg.*syp.*invOAF + sxg.*cyp;
    CondDelta = -cxg.*syp.*cos(zg-zp);
end

%% solution from YS

targetNed = [1;0;0];
dcmI2B = angle2dcm(0,-yr,-zr,'xyz');
rotmI2ShootNav = angle2dcm(0,0,-zp,'yxz');
targetNed = dcmI2B * targetNed;

z_hat = targetNed / norm(targetNed);
b_hat = [0,1,0];
b_hat = rotmI2ShootNav * b_hat';
y_hat = cross(z_hat,b_hat)/norm(cross(z_hat,b_hat));
x_hat = cross(y_hat,z_hat);
temp = [z_hat,x_hat,y_hat];
dcmNed2Image = temp';
[YS_Yaw,YS_Roll,YS_Pitch] = dcm2angle(dcmNed2Image,'zxy');

% rotmI2Gimbal = angle2dcm(-Pitch,-Roll,-Yaw,'yxz');
% gimbal_hat = [1,0,0];
% gimbal_result = rotmI2Gimbal * gimbal_hat';

%% Verification
testDcm = angle2dcm(-yg,-xg,-zg,'yxz');
testResult = zeros(3,n);
for i = 1:n
    testResult(:,i) = testDcm(:,:,i) * [1,0,0]';
end

%% Plotting
xg = xg*r2d; yg = yg*r2d; zg = zg*r2d;

figure(4)
clf
hold on
grid on
plot(zg,xg)
plot(zg,yg)
plot(zg,zg)
xlabel('Gimbal yaw angle')
ylabel('angle [deg]')
plot(YS_Yaw/d2r,YS_Pitch/d2r,'rx','MarkerSize',15,'LineWidth',5);
plot(YS_Yaw/d2r,YS_Roll/d2r,'rx','MarkerSize',15,'LineWidth',5);

yyaxis right
plot(zg,direcCond)
plot(zg,CondDelta);
title('directional condition','FontSize',15)
ylabel('directional condition value')
plot([-100 100],[0 0],'k')
legend('roll','pitch','yaw','directional condition','0 line','Location','northwest')
text(-20,0,0,'ㅋㅋㅋㅋㅋㅋㅋ ㅠㅠ','FontName','Times New Roman','fontsize',30)

figure(5)
clf
plot(sqrt(beta.^2./(1-alpha.^2)))

figure(6)
% verification plot
clf
grid on
hold on
plot3(testResult(1,:),testResult(2,:),testResult(3,:),'r.')
% plot3(testResult(2,:),testResult(1,:),-testResult(3,:),'ro')
% plot_arrow([0,0,0],[1,0,0],'k');
view(45,45)
axis equal
title('Did it 수렴?')

figure(7)
clf
grid on
hold on
plot(zg,direcCond)
title('directional condition','FontSize',15)