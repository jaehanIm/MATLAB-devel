clear all

% Simulation param.
simT = 1e2*3; % Simulation step
n = 5; % grid length
h = 1; % grid size
networkHeight = 3; %position of network
wallHeight = 1;
F = 0.5; % Simulation Speed

updateRate = 0.03;

% Dynamics param.
k_g = 1; % 사회적 인터넷 거부감 지표
k_n = 1;  % 사회적 인터넷 인프라 지표
c = .1;  % 사회 자율성 억제지표

agentIdx = 1:n^2;
agentPos = zeros(n^2,2); % agent position data
agentHeight = zeros(n^2,1); % agent vertical position
agentVel = zeros(n^2,1); % agent velocity position
agentSatisfaction = ones(n^2,1); % agent satisfaction
networkCenter = n*h/2*ones(2,1); % position of network
force = zeros(n^2,1); % force applied to each agent
networkConnectivity = 0; % How much are we connected to each other?
T = 1:simT;

history.pos = zeros(simT,n^2,1);
history.vel = zeros(simT,n^2,1);
history.satis = zeros(simT,n^2,1);
history.force = zeros(simT,n^2,1);

% agent pos. init.
temp = 1;
for i = 1:n
    for j = 1:n
        agentPos(temp,:) = [i*h-h/2,j*h-h/2];
        temp = temp+1;
    end
end

% agent desire init.
agentChar = rand(n^2,1); % nominal
% agentChar = rand(n^2,1)/10+0.9; % eager
% agentChar = rand(n^2,1)/30; % lousy

% agentSatisfaction = rand(n^2,1);
agentSatisfaction = agentSatisfaction * 0.4 + rand(n^2,1)*0.2;

prev_isWantingToConnect = zeros(n^2,1);

%% Main Simulation
for t=T
    for i = agentIdx
        % Update networkConnectivity
        networkConnectivity = networkConnectivityScore(agentHeight,networkHeight);
%         networkConnectivity = 1;
        % Update Agent Satisfaction
        if prev_isWantingToConnect(i) == 1
            agentSatisfaction(i) = agentSatisfaction(i) + (networkHeight-abs(agentHeight(i)-networkHeight)) * (1-agentChar(i)) * networkConnectivity * updateRate * F;
        else
            agentSatisfaction(i) = agentSatisfaction(i) - networkHeight * agentChar(i) * updateRate * F;
        end
        if agentSatisfaction(i) > 1
            agentSatisfaction(i) = 1;
        elseif agentSatisfaction(i) < 0
            agentSatisfaction(i) = 0;
        end
        % What do you want to do? Be alone? Want to be connected?
        isWantingToConnect = Roulette(agentSatisfaction(i));
        % Apply force
        force(i) = -k_g -c*agentVel(i);
        if isWantingToConnect
            force(i) = force(i) + (networkHeight - agentHeight(i))*k_n + k_g;
        end
        prev_isWantingToConnect(i) = isWantingToConnect;
    end
    for i = agentIdx
        % Update agent velocity
        agentVel(i) = agentVel(i) + force(i)*F;
        % Update agent position
        agentHeight(i) = agentHeight(i) + agentVel(i)*F;
        if agentHeight(i)<0
            agentHeight(i) = 0;
        elseif agentHeight(i) > networkHeight
            agentHeight(i) = networkHeight;
        end
        
    end
    history.pos(t,:,:) = agentHeight;
    history.vel(t,:,:) = agentVel;
    history.satis(t,:,:) = agentSatisfaction;
    history.force(t,:,:) = force;
    history.connectivity(t,:) = networkConnectivity;
end

%% Base plotting
figure(1)
clf
hold on
grid on

% Grid drawing
for i = 0:n
    plot3([i*h i*h],[0 n*h],[0 0],'k')
    plot3([0 n*h],[i*h i*h],[0 0],'k')
    plot3([i*h i*h],[0 n*h],[1 1],'k')
    plot3([0 n*h],[i*h i*h],[1 1],'k')
end

% Wall drawing
for i = 0:n
    for j = 0:n
        plot3([i*h i*h],[j*h j*h],[0 1],'k')
    end
end

% Network drawing
plot3(n*h/2,n*h/2,networkHeight,'ko','LineWidth',3)
plot3([n*h/2,n*h/2],[n*h/2,n*h/2],[0 networkHeight],'k:')

view(75,35)
axis equal
set(gcf, 'Position', [2100,100,500,500])

%% Data visualization
alpha = 0.5;
% figure(2)
% clf
% hold on
% for i = 1:n^2
%     plot(history.pos(:,i,1),'Color',[1,1,1]*alpha)
% %     plot(history.pos(:,i,1))
% end
% plot(mean(history.pos')+std(history.pos'),'b:','LineWidth',2)
% plot(mean(history.pos')-std(history.pos'),'b:','LineWidth',2)
% plot(mean(history.pos'),'r-','LineWidth',3)
% plot([0 simT],[networkHeight networkHeight],'k--','LineWidth',1)
% title('개인별 연결 위상')
% grid on
% ylim([0 3])
% 
% figure(3)
% clf
% hold on
% for i = 1:n^2
%    plot(1-history.satis(:,i,1),'Color',[1,1,1]*alpha)
% %    plot(1-history.satis(:,i,1))
% end
% plot(mean(1-history.satis')+std(1-history.satis'),'b:','LineWidth',2)
% plot(mean(1-history.satis')-std(1-history.satis'),'b:','LineWidth',2)
% plot(mean(1-history.satis'),'r-','LineWidth',3)
% title('연결 욕구')
% grid on
% ylim([0 1])
% 
% figure(4)
% clf
% hold on
% plot(history.connectivity(:)*100,'k')
% title('연결망 활성도')
% ylabel('%','FontSize',15)
% ylim([0 100])
% grid on

figure(5)
clf
hold on
subplot(3,1,1)
hold on
for i = 1:n^2
%     plot(history.pos(:,i,1)/3,'Color',[1,1,1]*alpha)
    plot(history.pos(:,i,1)/3)
end
plot(mean(history.pos'/3)+std(history.pos'/6),'b:','LineWidth',2)
plot(mean(history.pos'/3)-std(history.pos'/6),'b:','LineWidth',2)
plot(mean(history.pos'/3),'r-','LineWidth',3)
plot([0 simT],[networkHeight/3 networkHeight/3],'k--','LineWidth',1)
title('개인별 연결 위상')
grid on
ylim([0 networkHeight/3])

subplot(3,1,2)
hold on
for i = 1:n^2
%    plot(1-history.satis(:,i,1),'Color',[1,1,1]*alpha)
   plot(1-history.satis(:,i,1))
end
plot(mean(1-history.satis')+std(1-history.satis'),'b:','LineWidth',2)
plot(mean(1-history.satis')-std(1-history.satis'),'b:','LineWidth',2)
plot(mean(1-history.satis'),'r-','LineWidth',3)
title('연결 욕구')
grid on
ylim([0 1])

subplot(3,1,3)
hold on
plot(history.connectivity(:)*100,'k')
title('연결망 활성도')
ylabel('%','FontSize',15)
ylim([0 100])
grid on

%% Animation
figure(1)
aniConnectionLine = animatedline;
zlim([0 max(max(history.pos))])
blue = [0 0 1];
red = [1 0 0];
black = [0,0,0];
white = [1,1,1];

% World Plotting
for t = T(280:end)
    for i = 1:n^2
        % colored version
        aniAgent(i) = animatedline('Marker','o','MarkerSize',10,'LineStyle','none','Color',blue*(1-history.satis(t,i))+red*history.satis(t,i),'MarkerFaceColor',blue*(1-history.satis(t,i))+red*history.satis(t,i));
        % bw version
%         aniAgent(i) = animatedline('Marker','o','MarkerSize',10,'LineStyle','none');
        
        % colored line
        distanceFactor = abs(networkHeight - history.pos(t,i))/networkHeight;
%         aniLine(i) = animatedline('LineStyle','--','Color',distanceFactor * white + (1-distanceFactor) * red);
        addpoints(aniAgent(i),agentPos(i,1),agentPos(i,2),history.pos(t,i))
%         addpoints(aniLine(i),networkCenter(1),networkCenter(2),networkHeight)
%         addpoints(aniLine(i),agentPos(i,1),agentPos(i,2),history.pos(t,i))
    end
    drawnow
%     scene(t) = getframe(gcf);
    t
    pause(1)
    for i = 1:n^2
        clearpoints(aniAgent(i))
%         clearpoints(aniLine(i))
    end
end
% 
% v = VideoWriter('large.avi');
% v.FrameRate = 30;
% open(v);
% for i = 1:length(scene)
%     frame = scene(i);
%     writeVideo(v,frame);
% end
% 
% close(v);

