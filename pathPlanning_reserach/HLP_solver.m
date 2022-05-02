function [routeResult,score]=HLP_solver(map)
% require A C N ND vnum totLoad

% Parameter setting
V = map.vnum;
N = map.N-1; % number of superNodes
T = N + 1;
totLoad = map.totLoad;
ND = map.ND; % Super node demand matrix
ND_unit = ND;
for i = 1:T-1
    ND = horzcat(ND,ND_unit);
end

X = zeros(T,T,V); % selector
% Obg = map.obgSet; % obligatory edge pass constraint

% Map gen. / dist. calc.
Y = map.C;

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

loadConst = totLoad;
if V > 1 % vehicle demand equality constraint
    for v = 1:V
        X = zeros(T,T,V);
        X(:,:,v) = ND;
        A = vertcat(A,X(:)');
        loadConst = totLoad;
        b = vertcat(b,loadConst);
    end
end

% Equality Constraint
Aeq = [];
beq = [];

for j = 2:T % visiting constraint
    X = zeros(T,T,V);
    X(j,:,:) = 1;
    X(j,j,:) = 0;
    Aeq = vertcat(Aeq,X(:)');
    beq = vertcat(beq,1);
end

for v = 1:V % home depot constraint
    X = zeros(T,T,V);
    X(2:T,1,v) = 1;
%     X(1,2:T,v) = 1;
    Aeq = vertcat(Aeq,X(:)');
    beq = vertcat(beq,1);
end

for v = 1:V % inout constraint
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
Error = false;
flag = 1;
iterationNum = 0;

disp('initiating IP process!');
while flag == 1
    iterationNum = iterationNum + 1;
    options = optimoptions('intlinprog','AbsoluteGapTolerance',0.1,'IntegerTolerance',1e-6,'Display','none');
    [result,score,exitFlag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

    if exitFlag == -2
        disp("[HLP solver] ERROR : Unable to solve IP. No feasible solution existing");
        disp("ErrorVal : "+num2str(loadConst));
        Error = true;
        break;
    end

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
        tempResult = tempResult - diag(diag(tempResult)); %와 이거 찾느라 한세월 신발 
        pathResult = tempResult;
        break;
    end
    TEXT = ["Iteration Number : ",iterationNum];
    disp(TEXT);
end
disp('IP process complete!');
disp("loadConstVal : "+num2str(loadConst));

if Error
    routeResult = -1;
    score = -1;
else
    % Linker
    routeResult = detectRoutes(pathResult);
end

end