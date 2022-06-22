%% MC master for ACS work
clear all

addpath('./ACO')
addpath('./MILP')
addpath('./ComDetTBv090/')
addpath('./ComDetTBv090/Algorithms/')
addpath('./ComDetTBv090/Auxiliary/')
addpath('./ComDetTBv090/Graphs/')

vnum = 3;
depotPos = [60 60 60];
fovFactor = 3;
inpection_dist = 7; % Inspection distance
mapheight = 4.0;
conThres = 10;
antNo = 20;
stopThres = 200;
trialNums = 10;


MCData = [];
temp = [];
fovC = 1;
totC = 1;

fovFactorSet = [4:-0.5:1,0.75,0.5,0.4,0.3];
% conThresSet = 5:10:45;
conThresSet = [4,5,6,7,9,10,15,20,25,35,45];
% fovFactorSet = 1.5:-0.5:1.5;
% conThresSet = 40:5:40;
totalTestNum = length(fovFactorSet) * length(conThresSet);

for fovFactor = fovFactorSet
    conC = 1;
    for conThres = conThresSet
        completeTime_soleACO = [];
        disp("================")
        disp("case setting : ("+num2str(fovFactor)+", "+num2str(conThres)+")");
        disp("Test set progress : "+num2str(totC / totalTestNum * 100)+"% done");
        tryC = 1;
        disp("Generate map")
        poc_path_planner;
        disp("Initialization for ACS")
        MC_MAIN_soleACS;
        for tryNum = 1:trialNums
            clusteringTime = []; interCompleteTime = []; solveTime = []; totalScoreHistory = [];
            soleACO_time = []; soleACO_result = [];
            equiPerformanceTimeL = []; equiPerformanceTime = []; equiTimePerformance = [];
            disp("Try set progress : "+num2str(tryNum/trialNums*100)+"% of Test progress : "+num2str(totC / totalTestNum * 100)+"%")

            disp("MAIN start")
            MC_MAIN
            disp("soleACS start")
            ACSVRP_forSoleACS
    
            temp.TestClusteringTime = clusteringTime;
            temp.TestInterCompleteTime = interCompleteTime;
            temp.TestSolveTime = solveTime;
            temp.TestScoreHist = totalScoreHistory;
            temp.TestScoreHistL = totalScoreHistoryL;
            temp.TestSCoreHistPerV = vehScore;
    
            temp.CompCompleteTime = completeTime_soleACO;
            temp.CompSolveTime = soleACO_time;
            temp.CompScore = soleACO_result;
            temp.CompScoreL = soleACO_resultL;
            temp.CompScorePer = soleACO_resultPerV;
    
            temp.degreeConnectivity = degreeConnectivity;
            temp.N = N;

            temp.equiPerformanceTime = equiPerformanceTime;
            temp.equiPerformanceTimeL = equiPerformanceTimeL;
            temp.equiTimePerformance = equiTimePerformance;
    
            MCData{fovC, conC, tryC} = temp;
            tryC = tryC + 1;
        end
        conC = conC + 1;
        totC = totC + 1;
    end  
    fovC = fovC + 1;
end

MCData