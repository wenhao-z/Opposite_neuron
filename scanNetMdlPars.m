% Test the performance of decentralized system under different parameters
% Compared with scanNetMdlPars.m, this code considers network's estimate
% over time.

% Wen-Hao Zhang, Nov-17-2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

setWorkPath;

% Load parameters
% parsScanDoubleCANNs; % Two coupled CANNs
parsScanNetIntSeg; % Each net contains congruent (integration) and opposite (segregation) neurons
% parsDemoNetIntSeg;
% parsDemoNetIntSeg_diffAmpl;

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Net Simulation
NetStatStruct = struct('meanBumpPos', [], ...
    'mrlBumpPos', [], ...
    'concBumpPos', [], ...
    'varBumpPos', [], ...
    'OHeightAvg', []);
NetStatStruct = repmat(NetStatStruct, size(parGrid));

parpool(12);
tStart = clock;
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    if (netpars.Jrc + netpars.Jrp) > netpars.Jc
        % Parameters of persistent activity
        tmp = nan(netpars.numNets * netpars.numGroupPerNet, 1);
        for varName = fieldnames(NetStatStruct(iterPar))';
            NetStatStruct(iterPar).(varName{1}) = tmp;
        end      
    else
        % Network input     
        InputSet = makeNetInput([], netpars, struct('Iext', [], 'initRandStream', []));
        
        % Run simulation
        outArgs = struct('InputSet', [], ...
            'NetStat', NetStatStruct(iterPar));
        [InputSet, NetStatStruct(iterPar)] = simDecenNet(InputSet, netpars, outArgs);
    end
end

tEnd = clock;

%% Transform high-dim struct into a scalar struct with high-dim fields
for varName = fieldnames(NetStatStruct)';
    szFields = size(NetStatStruct(1).(varName{1}));
    NetStat = [NetStatStruct.(varName{1})];
    NetStat = reshape(NetStat, [szFields, size(parGrid)]);
    NetStat = permute(NetStat, [ndims(szFields)+1:ndims(NetStat), 1:ndims(szFields)]); % last dim is index of group
    NetEstimRes.(varName{1}) = squeeze(NetStat); % First several dims are the same as parGrids. The last few dims are the same as varName
end
clear varName NetStat NetStatStruct szFields

%% Bayes prediction
bayesRes = getBayesPrediction(NetEstimRes, NetPars, dimPar);

%% save
savePath = getSavePath(Path_RootDir, NetPars);
mkdir(savePath);

str = datestr(now, 'yymmddHHMM');
fileName = ['scanNetPars_', str(1:6), ...
    '_', str(7:end) '.mat'];

save(fullfile(savePath, fileName))
% save(fullfile(savePath, fileName), '-v7.3')