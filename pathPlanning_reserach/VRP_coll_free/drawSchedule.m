function [] = drawSchedule(figNum,N,testSchedule,vnum)

Color = rand(vnum,3);
Width = linspace(2,5,vnum);

figure(figNum)
clf
initTemp = zeros(N,N,N);
termTemp = zeros(N,N,N);
vehTemp = zeros(N,N,N);
checkEmpty = ones(N,N);
for i = 1:N
    for j = 1:N
        if ~isempty(testSchedule{i,j}.info)
            for k = 1:testSchedule{i,j}.num
                initTemp(i,j,k) = testSchedule{i,j}.info{k}(1);
                termTemp(i,j,k) = testSchedule{i,j}.info{k}(2);
                vehTemp(i,j,k) = testSchedule{i,j}.info{k}(3);
            end
            checkEmpty(i,j) = false;
        end
    end
end

for i = 1:N
    for j = 1:N
        if ~checkEmpty(i,j)
            for k = 1:testSchedule{i,j}.num
                plot([initTemp(i,j,k),termTemp(i,j,k)],[i,i],'LineWidth',Width(vehTemp(i,j,k)),'Color',Color(vehTemp(i,j,k),:))
                hold on
            end
        end
    end
end

end