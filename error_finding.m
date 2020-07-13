org = load('workspace.mat');
mod = load('mod_data_0412.mat');

close all

hold on
plot(org.time,org.posGP_ned(:,3)-mod.posGP_ned(:,3))
plot(org.time,org.posGP_ned(:,2)-mod.posGP_ned(:,2),'r')
plot(org.time,org.posGP_ned(:,1)-mod.posGP_ned(:,1),'k:')
figure()
hold on
plot(org.time,org.posGPS_ned(:,3)-mod.posGPS_ned(:,3))
plot(org.time,org.posGPS_ned(:,2)-mod.posGPS_ned(:,2),'r')
plot(org.time,org.posGPS_ned(:,1)-mod.posGPS_ned(:,1),'k:')