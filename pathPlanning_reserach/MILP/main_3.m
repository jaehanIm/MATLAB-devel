clear all

%% Optimization param setting

% Parameter setting
V = 3;
N = 15;
capacity = 20;
T = N + 1;

X = zeros(T,T,V); % selector

% Map gen. / dist. calc.
posX = 0;
posY = 0;
for i = 1:N
    posX(i+1) = rand(1) * 10 - 5;
    posY(i+1) = rand(1) * 10 - 5;
end

for i = 1:T
    for j = 1:T
        Y(i,j) = sqrt((posX(i)-posX(j))^2 + (posY(i)-posY(j))^2);
    end
end

% Cost function
funit = Y(:);
f = funit;
for i = 1:V-1
    f = vertcat(f,funit);
end

% Integer Constraint
intcon = 1:length(f);
intL = length(f);

% Inequality Constraint
A = [];
b = [];
for v = 1:V
    X = zeros(T,T,V);
    X(:,:,v) = 1;
    X = X.*Y;
    A = vertcat(A,X(:)');
    b = vertcat(b,capacity);
end

% Equality Constraint
Aeq = [];
beq = [];
for j = 2:T
    X = zeros(T,T,V);
    X(j,:,:) = 1;
    X(j,j,:) = 0;
    Aeq = vertcat(Aeq,X(:)');
    beq = vertcat(beq,1);
end

for v = 1:V
    X = zeros(T,T,V);
    X(2:T,1,v) = 1;
    Aeq = vertcat(Aeq,X(:)');
    beq = vertcat(beq,1);
end

for v = 1:V
    for j = 1:T
        X = zeros(T,T,V);
        X(j,1:T,v) = 1;
        X(1:T,j,v) = -1;
        X(j,j,v) = 0;
        Aeq = vertcat(Aeq,X(:)');
        beq = vertcat(beq,0);
    end
end


% Bound
lb = zeros(intL,1);
ub = ones(intL,1);

%% Run solver

% Solver
flag = 1;
iterationNum = 0;
while flag == 1
    iterationNum = iterationNum + 1;
    options = optimoptions('intlinprog','AbsoluteGapTolerance',0.1,'IntegerTolerance',0.001);
    result = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
    tempResult = reshape(result,[T,T,V]);
    tempResult = sum(tempResult,3);
    
    % Stagnation error compensation
    tempResult(tempResult<0.1) = 0;
    
    % Detect subtours
    subtour = detectSubtours(tempResult);
    flag = ~isempty(subtour); % check if there is any subtours generated
    
    % add constraint if required (subtour elimination constraint)
    if flag == 1
        for subNum = 1:size(subtour,1)
            curSubtour = subtour(subNum,:);
            curSubtour = curSubtour(find(curSubtour));
            X = zeros(T,T,V);
            X(curSubtour,curSubtour,1:V) = 1;
            for n = curSubtour
                X(n,n,:) = 0;
            end
            A = vertcat(A,X(:)');
            b = vertcat(b,size(curSubtour,2)-1);
        end
    else
        pathResult = tempResult;
        break;
    end
    TEXT = ["Iteration Number : ",iterationNum];
    disp(TEXT);
end

% Linker
routeResult = detectRoutes(pathResult);


%% Result plotting
figure(1)
clf
grid on
hold on
plot(posX(1),posY(1),'ro')
plot(posX(2:end),posY(2:end),'ko')
xlim([-5 5])
ylim([-5 5])

C = {'r','b','k','g','m'};
for v = 1:V
    pathData = routeResult(v,:);
    pathLen = length(find(pathData));
    for idx = 1:pathLen-1
        plot([posX(pathData(idx)) posX(pathData(idx+1))],[posY(pathData(idx)) posY(pathData(idx+1))],'color',C{v},'LineStyle','--');
    end
    plot([posX(pathData(idx+1)) posX(pathData(1))],[posY(pathData(idx+1)) posY(pathData(1))],'color',C{v},'LineStyle','--');
end