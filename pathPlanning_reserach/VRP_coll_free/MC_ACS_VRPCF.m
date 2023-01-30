
totalData = {};
vnumCount = 1;
distCount = 1;
MCiterNum = 100;
vnumCase = 4:12;
distThresCase = 5:10:45;
distThresCase = [5, 10, 13, 20, 25, 35];
% vnumCase = 12;
% distThresCase = 5;
tempRecord = zeros(100,1);
for distThres = distThresCase
    vnumCount = 1;
    for vnum = vnumCase
        dataRecord = [];
        for mm = 1:MCiterNum
            VRPCF_MAIN
            dataRecord(1,mm,:) = [VRPCF_solveTime,  VRPCF_result, graphDensity]; %보통 스코어 더 큼. 회피까지 해야하니까
            MAIN_sole_ACS
            dataRecord(2,mm,:) = [soleACO_time, soleACO_result, warn];
%             warn/soleACO_result/12
%             tempRecord(mm) = warn/soleACO_result/12;
        end
%         mean(tempRecord)
        totalData{distCount,vnumCount} = dataRecord;
        vnumCount = vnumCount+1;
    end
    distCount = distCount + 1;
end
save('totalData4.mat','totalData');
%% post analysis
load('totalData3.mat')
MCiterNum = 100;
vnumCase = 4:12;
distThresCase = [5, 10, 13, 20, 25, 35];
vCaseNum = length(vnumCase);
distCaseNum = length(distThresCase);

filteredTotalData = [];
for j = 1:vCaseNum
    dataRecord1 = [];
    dataRecord2 = [];
    dataRecord3 = [];
    dataRecord4 = [];
    dataRecord5 = [];
    for i = 1:distCaseNum
        localData = totalData{i,j};
        for k = 1:100
            if localData(1,k,3)<=3
                dataRecord1 = horzcat(dataRecord1,totalData{i,j}(:,k,:));
            elseif localData(1,k,3)<=5
                dataRecord2 = horzcat(dataRecord2,totalData{i,j}(:,k,:));
            elseif localData(1,k,3)<=10
                dataRecord3 = horzcat(dataRecord3,totalData{i,j}(:,k,:));
            elseif localData(1,k,3)<=15
                dataRecord4 = horzcat(dataRecord4,totalData{i,j}(:,k,:));
            else
                dataRecord5 = horzcat(dataRecord5,totalData{i,j}(:,k,:));
            end
        end
    end
    filteredTotalData{1,j} = dataRecord1;
    filteredTotalData{2,j} = dataRecord2;
    filteredTotalData{3,j} = dataRecord3;
    filteredTotalData{4,j} = dataRecord4;
    filteredTotalData{5,j} = dataRecord5;
end

totalData = filteredTotalData;

%% post analysis

vCaseNum = length(vnumCase);
distCaseNum = length(distThresCase);
filtered = zeros(vCaseNum, distCaseNum, 5); %vnum %[loss warnRatio]
% graphDensityFiltered = zeros(distCaseNum, vCaseNum, MCiterNum);
for i = 1:distCaseNum-1
    for j = 1:vCaseNum
        localData = totalData{i,j};
        localScoreLossRatio = mean((localData(1,:,2) - localData(2,:,2))./localData(2,:,2));
        localScoreLossRatioStd = std((localData(1,:,2) - localData(2,:,2))./localData(2,:,2));
        
        localWarnRatio = mean(localData(2,:,3)./localData(2,:,2)/vnumCase(j));
        localWarnRatioStd = std(localData(2,:,3)./localData(2,:,2)/vnumCase(j));
        
        lossToEarnRatio = localWarnRatio./localScoreLossRatio;
%         lossToEarnRatio(lossToEarnRatio<0) = 0;
        
        filtered(i,j,:) = [localScoreLossRatio, localWarnRatio, localScoreLossRatioStd, localWarnRatioStd, lossToEarnRatio];
%         graphDensityFiltered(i,j,:) = localData(1,:,3);
    end
end

markers = {'o','*','s','d','v','>','h'};

% filtered(2,:,:) = (filtered(2,:,:)+filtered(3,:,:))./2; % temporary
% filtered(3,:,:) = []; % temporary

figure(1)
clf
for i = 1:distCaseNum-1
    errorbar(vnumCase,filtered(i,:,1)*100,filtered(i,:,3)*10,'.--','Marker',markers{i})
    hold on
end
plot([4 12],[0 0],'k--')
title('Performance Gap','fontsize',13)
xlabel('Number of vehicles (V)','fontsize',13)
ylabel('Gap [%]','fontsize',13)
legend('$\mathcal{K}<3\%$','$\mathcal{K}<5\%$','$\mathcal{K}<10\%$','$\mathcal{K}<15\%$','$\mathcal{K}<20\%$','Location','best','Interpreter','latex','fontsize',12)
grid on
drawnow
exportgraphics(gca,'[0]Performance gap.png')

figure(2)
clf
for i = 1:distCaseNum-1
    errorbar(vnumCase,filtered(i,:,2)*100,filtered(i,:,4)*10,'.--','Marker',markers{i})
    hold on
end
title('Conflict Ratio','fontsize',13)
xlabel('Number of vehicles (V)','fontsize',13)
ylabel('Ratio [%]','fontsize',13)
legend('$\mathcal{K}<3\%$','$\mathcal{K}<5\%$','$\mathcal{K}<10\%$','$\mathcal{K}<15\%$','$\mathcal{K}<20\%$','Location','northwest','Interpreter','latex','fontsize',12)
grid on
drawnow
exportgraphics(gca,'[0]Conflict ratio.png')

% subplot(3,1,3)
% for i = 1:distCaseNum
%     plot(vnumCase,filtered(i,:,5),'.--')
%     hold on
% end

% figure(4)
% clf
% for i = 1:distCaseNum-1
%     temp = graphDensityFiltered(i,:,:);
%     temp = temp(:);
%     scatter(ones(size(temp))*i,temp)
%     hold on
% end

figure(3)
clf
plot(filtered(:,:,2)*100,filtered(:,:,1)*100,'k*')
grid on
xlabel('Conflict ratio [%]','fontsize',13)
ylabel('Performance gap [%]','fontsize',13)
title('Performance Gap to Conflict Ratio','fontsize',13)
drawnow
exportgraphics(gca,'[0]gap to ratio.png')

