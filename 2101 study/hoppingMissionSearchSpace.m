%% plotting search space
clear all
focalLengthSet = linspace(6.83,119.94,50);
InspectionDistanceSet = linspace(17,35,50);
time = [];
for focalLength = focalLengthSet
for inspectionDistance = InspectionDistanceSet
hoppingMissionDesign
% time = vertcat(time,[focalLength,inspectionDistance,reqShootTime]);
time = vertcat(time,[focalLength, inspectionDistance, reqTotalFlightTime]);
end
end
figure(1)
clf
grid on
hold on

for i = 1:length(time)
    if time(i,3) >= 20 * 60
        plot3(time(i,1),time(i,2),time(i,3)/60,'r.') 
    else
        plot3(time(i,1),time(i,2),time(i,3)/60,'k.')
    end
end
view(-45,45)
xlabel('focal Length [mm]')
ylabel('inspection distance [m]')
zlabel('inspection time [min]')