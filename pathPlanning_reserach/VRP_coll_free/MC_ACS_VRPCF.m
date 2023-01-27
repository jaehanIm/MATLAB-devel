
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

% filtered(distCaseNum,vCaseNum,test/comparison)
figure(1)
clf
subplot(2,1,1)
for i = 1:distCaseNum
    errorbar(vnumCase,filtered(i,:,1)*100,filtered(i,:,3)*10,'.-')
    hold on
end
plot([4 12],[0 0],'k--')
title('Performance degradation')
xlabel('Number of agents (N)')
ylabel('Degradation [%]')
legend('<5%','5~10%','10~20%','20~30%','30~ %','Location','best')
grid on

subplot(2,1,2)
for i = 1:distCaseNum
    errorbar(vnumCase,filtered(i,:,2)*100,filtered(i,:,4)*10,'.--')
    hold on
end
title('Conflict Ratio')
xlabel('Number of agents (N)')
ylabel('[%]')
legend('<5%','5~10%','10~20%','20~30%','30~ %','Location','best')
grid on

% subplot(3,1,3)
% for i = 1:distCaseNum
%     plot(vnumCase,filtered(i,:,5),'.--')
%     hold on
% end

figure(2)
clf
for i = 1:distCaseNum
    temp = graphDensityFiltered(i,:,:);
    temp = temp(:);
    scatter(ones(size(temp))*i,temp)
    hold on
end

figure(3)
clf
plot(filtered(:,:,2)*100,filtered(:,:,1)*100,'k*')
grid on
xlabel('conflict ratio [%]')
ylabel('Degradation [%]')
title('Performance degradation to conflict ratio')