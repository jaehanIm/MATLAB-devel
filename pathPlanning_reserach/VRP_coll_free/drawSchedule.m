function [] = drawSchedule(figNum,N,colony,vnum)

Color = rand(vnum,3);
Width = linspace(2,4,vnum);
testSchedule = colony.queen.reservation;
tour = colony.queen.tour;
tick = colony.queen.tickHistory;

figure(figNum)
clf

%% mode 1
% for i = 1:N
%     for j = 1:N
%         if ~isempty(testSchedule{i,j}.info)
%             for k = 1:testSchedule{i,j}.num
%                 initTemp(i,j,k) = testSchedule{i,j}.info{k}(1);
%                 termTemp(i,j,k) = testSchedule{i,j}.info{k}(2);
%                 vehTemp(i,j,k) = testSchedule{i,j}.info{k}(3);
%             end
%             checkEmpty(i,j) = false;
%         end
%     end
% end
% 
% for i = 1:N
%     for j = 1:N
%         if ~checkEmpty(i,j)
%             for k = 1:testSchedule{i,j}.num
%                 plot([initTemp(i,j,k),termTemp(i,j,k)],[i,i],'LineWidth',Width(vehTemp(i,j,k)),'Color',Color(vehTemp(i,j,k),:))
%                 hold on
%             end
%         end
%     end
% end

%% mode 2

for k = 1:10
    initTemp = zeros(N,N);
    termTemp = zeros(N,N);
    vehTemp = zeros(N,N);
    checkEmpty = ones(N,N);
    for i = 1:N
        for j = 1:N
            if ~isempty(testSchedule{i,j}.info) && testSchedule{i,j}.num >= k
                initTemp(i,j) = testSchedule{i,j}.info{k}(1);
                termTemp(i,j) = testSchedule{i,j}.info{k}(2);
                vehTemp(i,j) = testSchedule{i,j}.info{k}(3);
                checkEmpty(i,j) = false;
            end
        end
    end
    
    initTemp = initTemp(:);
    termTemp = termTemp(:);
    vehTemp = vehTemp(:);
    checkEmpty = checkEmpty(:);
    
    for i = 1:size(initTemp(:))
        if ~checkEmpty(i)
            plot([initTemp(i),termTemp(i)],[i,i],'LineWidth',Width(vehTemp(i)),'Color',Color(vehTemp(i),:))
            hold on
        end
    end
    grid on
end


% figure(figNum+1)
% clf
tickLen = size(tick,2);
for i = 1:tickLen-1
    for v = 1:vnum
        initIdx = tour(v,i);
        termIdx = tour(v,i+1);
        initTick = tick(v,i);
        termTick = tick(v,i+1);
        if initTick ~= 0 && termTick ~= 0
            linIdx = sub2ind([N,N],initIdx,termIdx);
            plot([initTick termTick],[linIdx,linIdx],'LineWidth',Width(v),'Color',Color(v,:),'LineStyle',':','Marker','diamond')
            hold on
        end
    end
end
% grid on

figure(figNum+1)
clf
tickLen = size(tick,2);
for i = 1:tickLen-1
    for v = 1:vnum
        initIdx = tour(v,i);
        termIdx = tour(v,i+1);
        initTick = tick(v,i);
        termTick = tick(v,i+1);
        if initTick ~= 0 && termTick ~= 0
            linIdx = sub2ind([N,N],initIdx,termIdx);
            plot([initTick termTick],[linIdx,linIdx],'LineWidth',Width(v),'Color',Color(v,:))
            hold on
        end
    end
end
grid on

end