% Read out integration probability from the joint activities of congruent
% and opposite neurons and then compare it with Bayesian prediction.

% Wen-Hao Zhang, Nov-09, 2017
% Carnegie Mellon University
% wenhaoz1@andrew.cmu.edu

setWorkPath;

% Particular parameters
NetPars.AmplRatio = 0.35*[1, 1]';
NetPars.AmplRatio = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet,1]);

% Cue 1 at 0 deg; cue 2 ranges from 0 to 180 deg.
NetPars.Posi = 0: 5: 180;
NetPars.Posi = [zeros(size(NetPars.Posi)); NetPars.Posi];
NetPars.Posi = repmat(NetPars.Posi, [NetPars.numNets, 1]);
NetPars.tLen = 200 * NetPars.tau;

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Simulation of Decision-making circuit
NetStatStruct = struct('OHeight', [], ...
    'BumpPos', [], ...
    'meanBumpPos', [], ...
    'mrlBumpPos', [], ...
    'concBumpPos', [], ...
    'varBumpPos', [], ...
    'OHeightAvg', []);
NetStatStruct = repmat(NetStatStruct, size(parGrid));
% InputSetArray = cell(size(parGrid));

% parpool(8);
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    % Network input
    InputSet = makeNetInput([], netpars);
    
    % Run simulation
    outArgs = struct('InputSet', [], ...
        'NetStat', NetStatStruct(iterPar));
    [InputSet, NetStatStruct(iterPar)] = simDecenNet(InputSet, netpars, outArgs);
    %     InputSetArray{iterPar} = InputSet.O;
end

%% Transform high-dim struct into a scalar struct with high-dim fields
for varName = fieldnames(NetStatStruct)';
    szFields = size(NetStatStruct(1).(varName{1}));
    NetStat = [NetStatStruct.(varName{1})];
    NetStat = reshape(NetStat, [szFields, size(parGrid)]);
    NetStat = permute(NetStat, [ndims(szFields)+1:ndims(NetStat), 1:ndims(szFields)]); % last dim is index of group
    NetEstimRes.(varName{1}) = squeeze(NetStat); % First several dims are the same as parGrids. The last few dims are the same as varName
end
clear varName NetStat NetStatStruct szFields

%% Read out integration probability from C&O neurons in combined cue condition
[ProbInt1, ProbInt_Bayes, alpha] = estimIntProb(NetEstimRes, NetPars, dimPar, 1);
ProbInt2 = estimIntProb(NetEstimRes, NetPars, dimPar, 2);
ProbInt3 = estimIntProb(NetEstimRes, NetPars, dimPar, 3);

% Figure
figure;

dX = diff(NetPars.Posi(1:2,:),1,1)';

subplot(2,2,1); hold on
plot(dX, NetEstimRes.OHeightAvg(:,1,1))
plot(dX, NetEstimRes.OHeightAvg(:,1,3))
legend('Congruent', 'Opposite')
xlabel('Cue disparity |x_1-x_2|')
ylabel('Firing rate (Hz)')
set(gca, 'xtick', 0:90:180)
axis tight; axis square
ylim([15, 30])

subplot(2,2,2); hold on
% cSpec = cool(size(NetPars.AmplRatio, 2));
plot(dX, squeeze(ProbInt1));
plot(dX, squeeze(ProbInt2));
plot(dX, squeeze(ProbInt3));
plot(dX, squeeze(ProbInt_Bayes), '--');
legend('Actual', 'Bayes')
xlabel('Cue disparity |x_1 - x_2|')
ylabel('Integration probability')
set(gca, 'xtick', 0:90:180, 'ytick', 0: 0.25:1)
axis tight; axis square
title(['alpha=', num2str(alpha)])


% Decision boundary with ratio of readout weight
subplot(2,2,3); hold on
wCong = 0.7: 0.01: 1.4; % The read out weight from C neurons, 
% the weight from O neuron is fixed at 1

xPoint = zeros(size(wCong));
for iter = 1: length(wCong)
    rateC = NetEstimRes.OHeightAvg(:,1,1);
    rateO = NetEstimRes.OHeightAvg(:,1,3);
    
    funHandle = @(x) abs(interp1(dX, rateC*wCong(iter), x) -...
        interp1(dX, rateO, x));
    
    [~, Idx] = min(abs(rateC*wCong(iter) - rateO));
    xPoint(iter) = fminsearch(funHandle, dX(Idx));
end
plot(wCong, xPoint)
xlabel('W_{cong}/W_{oppo}')
ylabel('Decision boundary')
axis tight; axis square
set(gca, 'ytick', 0:90:180, 'xtick', 0.7:0.1: 1.4)

plot([wCong(1), 1], 90*ones(1,2), '--k')
plot(ones(1,2), [0, 90], '--k')