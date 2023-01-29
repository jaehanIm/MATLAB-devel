
totalData = {};
vnumCount = 1;
distCount = 1;
MCiterNum = 100;
vnumCase = 4:12;
distThresCase = 5:10:45;
distThresCase = [5, 10, 15, 20, 25, 35];
for distThres = distThresCase
    vnumCount = 1;
    for vnum = vnumCase
        dataRecord = [];
        for mm = 1:MCiterNum
            VRPCF_MAIN
            dataRecord(1,mm,:) = [VRPCF_solveTime,  VRPCF_result, graphDensity]; %보통 스코어 더 큼. 회피까지 해야하니까
            MAIN_sole_ACS
            dataRecord(2,mm,:) = [soleACO_time, soleACO_result, warn];
        end
        totalData{distCount,vnumCount} = dataRecord;
        vnumCount = vnumCount+1;
    end
    distCount = distCount + 1;
end
%% post analysis

vCaseNum = length(vnumCase);
distCaseNum = length(distThresCase);
filtered = zeros(vCaseNum, distCaseNum, 5); %vnum %[loss warnRatio]
graphDensityFiltered = zeros(distCaseNum, vCaseNum, MCiterNum);
for i = 1:distCaseNum
    for j = 1:vCaseNum
        localData = totalData{i,j};
        localScoreLossRatio = mean((localData(1,:,2) - localData(2,:,2))./localData(2,:,2));
        localScoreLossRatioStd = std((localData(1,:,2) - localData(2,:,2))./localData(2,:,2));
        
        localWarnRatio = mean(localData(2,:,3)./localData(2,:,2)/vnumCase(i));
        localWarnRatioStd = std(localData(2,:,3)./localData(2,:,2)/vnumCase(i));
        
        lossToEarnRatio = localWarnRatio./localScoreLossRatio;
%         lossToEarnRatio(lossToEarnRatio<0) = 0;
        
        filtered(i,j,:) = [localScoreLossRatio, localWarnRatio, localScoreLossRatioStd, localWarnRatioStd, lossToEarnRatio];
        graphDensityFiltered(i,j,:) = localData(1,:,3);
    end
end

filtered(2,:,:) = (filtered(2,:,:)+filtered(3,:,:))./2;
filtered(3,:,:) = [];
figure(1)
clf
for i = 1:distCaseNum-1
    errorbar(vnumCase,filtered(i,:,1)*100,filtered(i,:,3)*10,'.-')
    hold on
end
plot([4 12],[0 0],'k--')
title('Performance Gap','fontsize',15)
xlabel('Number of agents (N)','fontsize',15)
ylabel('Gap [%]','fontsize',15)
legend('<5%','<10%','<15%','<20%','>20%','Location','best','fontsize',13)
grid on
drawnow
exportgraphics(gca,'[0]Performance gap.png')

figure(2)
clf
for i = 1:distCaseNum-1
    errorbar(vnumCase,filtered(i,:,2)*100,filtered(i,:,4)*10,'.--')
    hold on
end
title('Conflict Ratio','fontsize',15)
xlabel('Number of agents (N)','fontsize',15)
ylabel('[%]','fontsize',15)
legend('<5%','<10%','<15%','<20%','>20%','Location','northwest','fontsize',13)
grid on
drawnow
exportgraphics(gca,'[0]Conflict ratio.png')

% subplot(3,1,3)
% for i = 1:distCaseNum
%     plot(vnumCase,filtered(i,:,5),'.--')
%     hold on
% end

figure(4)
clf
for i = 1:distCaseNum-1
    temp = graphDensityFiltered(i,:,:);
    temp = temp(:);
    scatter(ones(size(temp))*i,temp)
    hold on
end

figure(3)
clf
plot(filtered(:,:,2)*100,filtered(:,:,1)*100,'k*')
grid on
xlabel('Conflict ratio [%]','fontsize',15)
ylabel('Gap [%]','fontsize',15)
title('Performance Gap to Conflict Ratio','fontsize',15)
drawnow
exportgraphics(gca,'[0]gap to ratio.png')