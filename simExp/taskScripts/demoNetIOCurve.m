% Plot the input-firing rate curves of network model
% Each network module contains both congruent and opposite neurons.

% Author: Wen-Hao Zhang, June-12, 2018
% wenhao.zhang@pitt.edu
% University of Pittsburgh

% Particular parameters
NetPars.AmplRatio = 0:0.05:1.5;

NetPars.cueCond = 0;
NetPars.tLen = 20 * NetPars.tau;
NetPars.Posi = 0;
NetPars.Posi = repmat(NetPars.Posi, [NetPars.numNets*NetPars.numGroupPerNet, 1]);


% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);
%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Run simulation
NetStatStruct = struct('OAvgXTime', [], ...
    'OStdXTime', [], ...
    'OHeightAvg', []);
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
figure; hold on;
cSpec = [42, 194, 64;
    44, 40, 187;
    255, 138, 55;]/255;
cSpecOppo = 1 - cSpec;

IdxNeuron = 45; %90;

for iterGroup = [1,3]
    if iterGroup <= 2
        % Congruent cell
        colorLine = cSpec;
    else
        % Opposite cell
        colorLine = cSpecOppo;
    end
    
    plot(NetPars.AmplRatio, NetEstimRes.OHeightAvg(:, iterGroup))
end

axis tight; 
set(gca, 'xtick', 0:0.3:1.5, 'ytick', 0:10:50, 'fontsize', 12, 'ylim', [0, 55])
xlabel('Itensity of Cue (U_0)'); 
ylabel('Firing rate (Hz)')

legend('Congruent', 'Opposite', 'location', 'best')
title(['JrcRatio=', num2str(NetPars.JrcRatio), ...
', JrpRatio=', num2str(NetPars.JrpRatio), ...
', krpRatio=', num2str(NetPars.krpRatio)])
