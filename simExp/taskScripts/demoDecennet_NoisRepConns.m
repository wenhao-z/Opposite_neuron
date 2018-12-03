% Demo of the decentralized network model with congruent and opposite
% neurons with the reciprocal connections being added with noises to
% produce similar distributions of cong. and oppo. neurons as indicated in
% Gu et al., JNS 2006.

% To test the difference of tuning preference under each single cues, we
% need to scan the tuning curves of every neurons.

% Wen-Hao Zhang, July-5, 2017
% wenhaoz1@andrew.cmu.edu
% Carnegie Mellon University

%% Particular parameters for this task
NetPars.AmplRatio = [1, 1]';
NetPars.AmplRatio = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet,1]);

NetPars.cueCond = 1:2;
NetPars.tLen = 20 * NetPars.tau;
NetPars.Posi = -180: 10: 180;
NetPars.Posi(1) = [];
NetPars.Posi = repmat(NetPars.Posi, [NetPars.numNets*NetPars.numGroupPerNet, 1]);

datPath = fullfile(Path_RootDir, 'Data');
Folder = 'NetIntSeg';
fileName = 'scanNetPars_170705_1819.mat';
NetParsScan = load(fullfile(datPath, Folder, fileName), 'NetPars');

NetPars.stdWrepNoisRatio = 0.9; % Relative to the maximal weight of reciprocal connections
% NetPars.seedWrepNois = sum(clock) * 100;
NetPars.seedWrepNois = NetParsScan.NetPars.seedWrepNois;

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Run simulation
NetStatStruct = struct('OAvgXTime', [], ...
    'OStdXTime', []);
NetStatStruct = repmat(NetStatStruct, size(parGrid));

tStart = clock;
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/%d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    % Network input
    InputSet = makeNetInput([], netpars);
    %     InputSet.Iext = cat(3, InputSet.Iext, zeros(size(InputSet.Iext)));
    
    % Run simulation
    outArgs = struct('InputSet', [], ...
        'NetStat', NetStatStruct(iterPar));
    [InputSet, NetStatStruct(iterPar)] = simDecenNet(InputSet, netpars, outArgs);
end
tEnd = clock;

%% Transform high-dim struct into a scalar struct with high-dim fields
for varName = fieldnames(NetStatStruct)';
    szFields = size(NetStatStruct(1).(varName{1}));
    if (strcmp(varName, 'BumpPos') || strcmp(varName, 'OHeight')) && NetPars.nTrials >1
        NetStat = {NetStatStruct.(varName{1})};
        NetStat = cellfun(@(x) reshape(x(:,NetPars.tStat/NetPars.dt+1:end,:), size(x,1), []), NetStat, 'uniformout', 0);
        szFields = size(NetStat{1});
        NetStat = cell2mat(NetStat);
    else
        NetStat = [NetStatStruct.(varName{1})];
    end        
    NetStat = reshape(NetStat, [szFields, size(parGrid)]);
    NetStat = permute(NetStat, [ndims(szFields)+1:ndims(NetStat), 1:ndims(szFields)]); % last dim is index of group
    NetEstimRes.(varName{1}) = NetStat;
end
clear varName NetStat NetStatStruct szFields

%% Estimate the preferred directions under two single cues

cirPos = exp(1i * NetPars.Posi(1,:)*pi/NetPars.Width)';
tuneCurvePos = mean(bsxfun(@times, cirPos, NetEstimRes.OAvgXTime), 1);
tuneCurvePos = angle(tuneCurvePos)*NetPars.Width/pi;
tuneCurvePos = squeeze(tuneCurvePos); % [cueCond, Neuron index, Group index]

diffCurvePos_singleCues = abs(diff(tuneCurvePos(1:2,:,:), 1, 1));
Idx = (diffCurvePos_singleCues>NetPars.Width);
diffCurvePos_singleCues(Idx) = 2*NetPars.Width - diffCurvePos_singleCues(Idx);
diffCurvePos_singleCues = squeeze(diffCurvePos_singleCues);
clear Idx

%% Plot 
% Edge = linspace(0, NetPars.Width, 10);
% histDiffCurvePos = histc(diffCurvePos_singleCues(:), Edge);
% bar(Edge, histDiffCurvePos);

nEdge = 10;
edgeDiff = linspace(0, 180, nEdge + 1);
IdxGroup = [1,3; 2,4];
for iter = 1: 2
    subplot(1,2,iter)
    hist(reshape(diffCurvePos_singleCues(:,IdxGroup(iter,:)),1,[]), 10)
    % histDiffNet1 = histcounts(reshape(diffCurvePos_singleCues(:,[1,3]),1,[]), edgeDiff);
    % bar()
    xlabel('Difference of preferred direction')
    ylabel('Number of neurons')
    axis square; axis tight;
    set(gca, 'xlim', [-2, 180+2], 'xtick', 0:45:180, 'ylim', [0, 120], 'ytick', 0:25:100)
    title(['Module',  num2str(iter)])
end

% subplot(1,2,2)
% hist(reshape(diffCurvePos_singleCues(:,[2,4]),1,[]), 10)
% xlabel('Difference of preferred direction')
% ylabel('Number of neurons')
% axis square; axis tight;
% set(gca, 'xlim', [0, 180], 'xtick', 0:45:180, 'ylim', [0, 120], 'ytick', 0:25:100)
% title('Module 2')

