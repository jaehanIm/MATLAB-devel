function [] = drawBestTour_forSoleACO(colony , graph)

global vehNum homeIdx

figure(1)
hold on
queenTour = colony.queen.tour;
queenTourLen = colony.queen.vehTourLen;
hold on
%
for i = 1:vehNum
    color = rand(1,3);
    for j = 1:size(queenTour(i,1:queenTourLen(i)),2)  
        currentNode = queenTour(i,j);
        nextNode = queenTour(i,j+1);
        
        x1 = graph.node.x(currentNode);
        y1 = graph.node.y(currentNode);
        z1 = graph.node.z(currentNode);

        x2 = graph.node.x(nextNode);
        y2 = graph.node.y(nextNode);
        z2 = graph.node.z(nextNode);

        X = [x1 , x2];
        Y = [y1, y2];
        Z = [z1, z2];
        plot3(X, Y, Z+1,'Color',color,'LineWidth',3);
    end
end

plot3(graph.node.x(homeIdx),graph.node.y(homeIdx),graph.node.z(homeIdx)+1,'ok', 'markerSize' , 10 , 'MarkerEdgeColor' , 'r' , 'MarkerFaceColor', 'g');
for i = 2 : graph.n
    
    X = [graph.node.x(i)];
    Y = [graph.node.y(i)];
    Z = [graph.node.z(i)];
    
    plot3(X, Y, Z, 'ok', 'markerSize' , 3 , 'MarkerEdgeColor' , 'r' , 'MarkerFaceColor', [1, 0.6, 0.6]);
end

title('Route Result')
box('on');
grid on
% xlim([0 28.1])
% ylim([0 24.6])
% ylim([-inf inf])

%



end