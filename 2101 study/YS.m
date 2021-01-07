d2r = pi/180;
targetNed = [25,0,0];
dcmI2B = angle2dcm(0,30*d2r,-10*d2r,'xyz');
targetNed = dcmI2B * targetNed';
z_hat = targetNed / norm(targetNed);
b_hat = [0,1,0];
y_hat = cross(z_hat,b_hat)/norm(cross(z_hat,b_hat));
x_hat = cross(y_hat,z_hat);
x_hat = x_hat'; y_hat = y_hat';
temp = [z_hat,x_hat,y_hat];
dcmNed2Image = temp';

[Yaw,Roll,Pitch] = dcm2angle(dcmNed2Image,'zxy');
[Pitch,Roll,Yaw]/d2r
rotmI2Gimbal = angle2dcm(-Pitch,-Roll,-Yaw,'yxz');
gimbal_hat = [1,0,0];
gimbal_result = rotmI2Gimbal * gimbal_hat';

gimbal_result
z_hat