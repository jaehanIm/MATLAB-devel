%script to add all Community Detection Toolbox folders to MATLAB path
a = mfilename('fullpath');
a = a(1:end-8);
addpath(a)
addpath([a '/Algorithms'])
addpath([a '/Auxiliary'])
addpath([a '/Cluster_Number'])
addpath([a '/ClusterValidity'])
addpath([a '/Evaluation'])
addpath([a '/Experiments'])
addpath([a '/Graphs'])
addpath([a '/Help'])
%savepath