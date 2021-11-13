function [ colony ] = createColony( graph, colony , antNo, tau, eta, alpha,  beta)

global exploreP
nodeNo = graph.n;

for i = 1 : antNo
    
    initial_node = randi( [1 , nodeNo] ); % select a random node 
    colony.ant(i).tour(1) = initial_node;
    
    for j = 2 : nodeNo % to choose the rest of nodes 
        if rouletteWheel([exploreP 1-exploreP]) == 2
            
           currentNode =  colony.ant(i).tour(end);

           P_allNodes = tau( currentNode , :  ) .^ alpha .* eta( currentNode , :  )  .^ beta;
           P_allNodes(colony.ant(i).tour) = 0 ;  % Passing 0 to all the nodes visited so far
           P = P_allNodes ./ sum(P_allNodes);

           nextNode = rouletteWheel(P); 
           colony.ant(i).tour = [  colony.ant(i).tour , nextNode ];
           
        else
            
            currentNode =  colony.ant(i).tour(end);
            
            P_allNodes = eta( currentNode , :  )  .^ beta;
%             P_allNodes = ones(size(tau));
            P_allNodes(colony.ant(i).tour) = 0 ;  % Passing 0 to all the nodes visited so far
            P = P_allNodes ./ sum(P_allNodes);

            nextNode = rouletteWheel2(P); 
            colony.ant(i).tour = [  colony.ant(i).tour , nextNode ];
            
        end
    end
    
    % complete the tour 
    colony.ant(i).tour = [ colony.ant(i).tour , colony.ant(i).tour(1)]; %back to home!
end

end