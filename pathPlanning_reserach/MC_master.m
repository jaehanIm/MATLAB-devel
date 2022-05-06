%% MC master for ACS work

clear all

addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')
addpath('./ComDetTBv090/Graphs/')

vnum = 4;
depotPos = [60 60 60];
fovFactor = 3;
inpection_dist = 7; % Inspection distance
mapheight = 4.0;
conThres = 10;
antNo = 20;
stopThres = 400;
trialNums = 2;


MCData = [];
temp = [];
fovC = 1;
totC = 1;

fovFactorSet = 3:-0.5:0.5;
conThresSet = 5:5:30;
totalTestNum = length(fovFactorSet) * length(conThresSet);

for fovFactor = fovFactorSet
    conC = 1;
    for conThres = conThresSet
        disp("================")
        disp("case setting : ("+num2str(fovFactor)+", "+num2str(conThres)+")");
        disp("Test set progress : "+num2str(totC / totalTestNum * 100)+"% done");
        tryC = 1;
        for tryNum = 1:trialNums
            clusteringTime = []; interCompleteTime = []; solveTime = []; totalScoreHistory = [];
            completeTime_soleACO = []; soleACO_time = []; soleACO_result = [];
            disp("Try set progress : "+num2str(tryNum/trialNums*100)+"% of Test progress : "+num2str(totC / totalTestNum * 100)+"%")
            disp("MAIN start")
            MC_MAIN
            disp("soleACS start")
            MC_MAIN_soleACS
    
            temp.TestClusteringTime = clusteringTime;
            temp.TestInterCompleteTime = interCompleteTime;
            temp.TestSolveTime = solveTime;
            temp.TestScoreHist = totalScoreHistory;
    
            temp.CompCompleteTime = completeTime_soleACO;
            temp.CompSolveTime = soleACO_time;
            temp.CompScore = soleACO_result;
    
            temp.degreeConnectivity = degreeConnectivity;
            temp.N = N;
    
            MCData{fovC, conC, tryC} = temp;
            tryC = tryC + 1;
        end
        conC = conC + 1;
        totC = totC + 1;
    end  
    fovC = fovC + 1;
end

MCData