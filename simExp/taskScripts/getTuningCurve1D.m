% Get the 1D tuning curve of neurons under three stimulus conditions
% Decentralized neural network with two newtork modules
% Each network module contains both congruent and opposite neurons.

% Author: Wen-Hao Zhang, June-2, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Particular parameters
NetPars.AmplRatio = [0.35, 0.8]';
% NetPars.AmplRatio = [1, 1]';
NetPars.AmplRatio = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet,1]);

% NetPars.cueCond = 1:2;
NetPars.tLen = 20 * NetPars.tau;
NetPars.Posi = -180: 15: 180;
NetPars.Posi = repmat(NetPars.Posi, [NetPars.numNets*NetPars.numGroupPerNet, 1]);

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
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
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
    NetStat = [NetStatStruct.(varName{1})];
    NetStat = reshape(NetStat, [szFields, size(parGrid)]);
    NetStat = permute(NetStat, [ndims(szFields)+1:ndims(NetStat), 1:ndims(szFields)]); % last dim is index of group
    NetEstimRes.(varName{1}) = NetStat;
end
clear varName NetStat NetStatStruct szFields


%% Plot
figure
cSpec = [42, 194, 64;
    44, 40, 187;
    255, 138, 55;]/255;
cSpecOppo = 1 - cSpec;

IdxNeuron = 45; %90;
textTitle = {'Net 1, congruent', 'Net 2, congruent', ...
    'Net 1, opposite','Net 2, opposite'};
hAxe = [];

for iterGroup = 1: 4
    hAxe = [hAxe, subplot(2, 2, iterGroup)];
    hold on
    if iterGroup <= 2
        % Congruent cell
        colorLine = cSpec;        
    else
        % Opposite cell
        colorLine = cSpecOppo;
    end
    for cueCond = 1: 3
        hPlot(cueCond) = plot(NetPars.Posi(1,:), NetEstimRes.OAvgXTime(:, cueCond, IdxNeuron, iterGroup), ...
            '-', 'color', colorLine(cueCond,:), 'linew', 1.5);
                errorbar(NetPars.Posi(1,:), NetEstimRes.OAvgXTime(:, cueCond, IdxNeuron, iterGroup), ...
                    NetEstimRes.OStdXTime(:, cueCond, IdxNeuron, iterGroup), ...
                    'color', colorLine(cueCond,:), 'linew', 1.5)
%         fill([NetPars.Posi(1,:), flip(NetPars.Posi(1,:), 2)], ...
%             [NetEstimRes.OAvgXTime(:, cueCond, IdxNeuron, iterGroup)' + NetEstimRes.OStdXTime(:, cueCond, IdxNeuron, iterGroup)', ...
%             flip(NetEstimRes.OAvgXTime(:, cueCond, IdxNeuron, iterGroup)' - NetEstimRes.OStdXTime(:, cueCond, IdxNeuron, iterGroup)', 2)], ...
%             lines(1), 'linestyle', 'none', 'facecolor', colorLine(cueCond,:), ...
%             'facealpha', 0.3)
    end
    axis tight; axis square
    set(gca, 'xtick', -180: 90: 180, 'fontsize', 12, 'ylim', [0, 55])
    xlabel('Cue'); ylabel('Firing rate')
    title(textTitle{iterGroup})
    legend(hPlot, 'Combined', 'Cue 1 ', 'Cue 2')
end

axes(hAxe(1))
text(180, 20, ['\alpha_1=' num2str(NetPars.AmplRatio(1)) ',\alpha_2=' num2str(NetPars.AmplRatio(2))], ...
    'horizontalalignment', 'right')
