% Simulate the experiments to generate the results of decentralized network
% model consisting of congruent and opposite neurons.

% Author: Wen-Hao Zhang, June-2, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load parameters
setWorkPath;
addpath(fullfile(Path_RootDir, 'simExp'));
addpath(fullfile(Path_RootDir, 'simExp', 'taskScripts'));

parsDemoNetIntSeg; % Each net contains congruent (integration) and opposite (segregation) neurons

%% Perform virtual experiments
flagTask = 3;
% 1. 1D tuning curves under three cueing conditions
% 2. 2D tuning curves under three cueing conditions
% 3. Population activities under three cueing conditions,
%    and the decoded stimulus (bump position) under three cueing
%    conditions
% 4. Discrimination task of whether x_1 is larger than 0 deg.
% 5. Discrimination task of whether x_1 is larger than x_2.
% 6. Adding noises onto reciprocal connections across network modules to
%    produce similar distribution of cong. and oppo. neurons as indicated
%    in Gu et al., JNS 2006.
% 7. Read out integration probability from C and O neurons and compare it
%    with Bayesian prediction
% 8. Plot the input-firing rate curves of network model

switch flagTask
    case 1
        getTuningCurve1D;
    case 2
        getTuningCurve2D;
    case 3 
        demoDecenNet_PopRes;
    case 4
        typeDiscrimTask = 1;
        getNeuroMetricFunc;
    case 5
        typeDiscrimTask = 2;
%         getNeuroMetricFunc;
        newDiscrimTask;
    case 6
        demoDecennet_NoisRepConns;
    case 7
        demoReadIntProb;
    case 8
        demoNetIOCurve;
end
