function [ ] = drawBestTour2(colony , graph)
figure(11)
clf
queenTour = colony.queen.tour;
hold on
for i = 1 : length(queenTour) - 1
    
    currentNode = queenTour(i);
    nextNode =  queenTour(i+1);
    
    x1 = graph.node.x(currentNode);
    y1 = graph.node.y(currentNode);
    z1 = graph.node.z(currentNode);
    
    x2 = graph.node.x(nextNode);
    y2 = graph.node.y(nextNode);
    z2 = graph.node.z(nextNode);
    
    X = [x1 , x2];
    Y = [y1, y2];
    Z = [z1, z2];
    plot3(X, Y, Z, '-r');

end


plot3(graph.node.x(1),graph.node.y(1),graph.node.z(1),'ok', 'markerSize' , 10 , 'MarkerEdgeColor' , 'r' , 'MarkerFaceColor', 'g');

for i = 2 : graph.n
    
    X = [graph.node.x(i)];
    Y = [graph.node.y(i)];
    Z = [graph.node.z(i)];
    
    plot3(X, Y, Z, 'ok', 'markerSize' , 10 , 'MarkerEdgeColor' , 'r' , 'MarkerFaceColor', [1, 0.6, 0.6]);
end

title('Best tour (the queen)')
box('on');
grid on