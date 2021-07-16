function [ ] = drawBestTour2(colony , graph)
figure(11)
clf
queenTour = colony.queen.tour;
hold on
for i = 1 : length(queenTour) - 1
    
    currentNode = queenTour(i)
    nextNode =  queenTour(i+1);
    
    x1 = graph.node.x(currentNode)
    y1 = graph.node.y(currentNode);
    
    x2 = graph.node.x(nextNode);
    y2 = graph.node.y(nextNode);
    
    X = [x1 , x2];
    Y = [y1, y2];
    plot (X, Y, '-r');

end


for i = 1 : graph.n
    
    X = [graph.node.x(:)];
    Y = [graph.node.y(:)];
    
    plot(X, Y, 'ok', 'markerSize' , 10 , 'MarkerEdgeColor' , 'r' , 'MarkerFaceColor', [1, 0.6, 0.6]);
end

title('Best tour (the queen)')
box('on');