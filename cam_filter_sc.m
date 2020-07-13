time_sim = seconds(time);

velcam_sim_0 = timeseries(data.gimbalRPY_0,time_sim);
velcam_sim_1 = timeseries(data.gimbalRPY_1,time_sim);
cam_filter
close all
figure()
plot(out.test_u)
hold on
grid on
plot(data.gimbalRPY_1)
plot(velNav(:,1))