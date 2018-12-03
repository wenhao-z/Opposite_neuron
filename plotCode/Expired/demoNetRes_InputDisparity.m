% Demo of the responses of decentralized system 
%  for info. integration and segregation.
% The decentralized system is composed of several reciprocally connected
% network modules;
% Each network contains congruent and opposite neurons.

% Author: Wen-Hao Zhang, 5-April-2016

bPlotIntRes = 1;
flagStimDim = 1; % Dim of stimulus features

% set the parameters of the network
parsNetMdl2;

% ------- Adjustable parameters ---------------
NetPars.Ampl = 1.2*ones(2,1) * NetPars.U0;
NetPars.Jrc = 0.5*NetPars.Jc;
NetPars.JrpRatio = 0.7;
NetPars.cueCond = 0;
NetPars.tLen = 50 * NetPars.tau;

switch flagStimDim
    case 1
        Posi = -180: 20: 180;
%         Posi = 0: 45: 180;
        Posi = [zeros(size(Posi)); Posi];
%         Posi = flipud(Posi);
    case 2
        Posi = -180: 45: 180;
        [Posi1, Posi2] = meshgrid(Posi, Posi);
        Posi = [Posi1(:), Posi2(:)]';
end
% Posi = [0;180];
NetPars.Posi = repmat(Posi, [NetPars.numNets, 1]);
clear Posi Posi1 Posi2

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

%% Run simulation
OAvgArray = zeros(NetPars.N, 2*NetPars.numNets, numel(parGrid)); % [1c, ..., numNets c, 1o,..., numNets o]
meanEstimNetSim = cell(size(parGrid)); 
varEstimNetSim = cell(size(parGrid));
concEstimNetSim = cell(size(parGrid));
meanEstimComb = cell(size(parGrid));
varEstimComb = cell(size(parGrid));

for iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    % Load parameters
    netpars = parGrid(iterPar);
    
    % Network input
    netpars.Ampl = repmat(netpars.Ampl, [netpars.numNets, 1]);
    InputSet = makeNetInput(netpars.Posi, netpars);
    
    % Run simulation
    InputSet = simDecenNet2(InputSet, netpars);
    OAvg = mean(InputSet.O(:,:,1001:end), 3); % Temporal average
    OAvgArray(:,:,iterPar) = OAvg;
    
    % A linear read out 
    InputSet.OReadOut = sum(InputSet.O(:,1:2:end,:) , 2);
    InputSet.OReadOut(:,2,:) = sum(InputSet.O(:, 2:2:end,:) , 2);

    % Estimate bump position by population vector
    [BumpPos, meanBumpPos, varBumpPos, concBumpPos]  ...
        = getBumpPos(InputSet.O(:,:,1001:end), NetPars);
    meanEstimNetSim{iterPar} = shiftdim(meanBumpPos, -ndims(parGrid)+1);
    varEstimNetSim{iterPar} = shiftdim(varBumpPos, -ndims(parGrid)+1);
    concEstimNetSim{iterPar} = shiftdim(concBumpPos, -ndims(parGrid)+1);
    
    [~, meanBumpPos, varBumpPos]  = getBumpPos(InputSet.OReadOut(:,:,1001:end), NetPars);
    meanEstimComb{iterPar} = shiftdim(meanBumpPos, -ndims(parGrid)+1);
    varEstimComb{iterPar} = shiftdim(varBumpPos, -ndims(parGrid)+1);
    
    clear meanBumpPos varBumpPos
    
    if bPlotIntRes
        %% Plot results
        plotIntResDemo;        
        pause
    end
end

meanEstimNetSim = squeeze(cell2mat(meanEstimNetSim)); % the last dimension is index of network
varEstimNetSim = squeeze(cell2mat(varEstimNetSim));
meanEstimComb = squeeze(cell2mat(meanEstimComb)); % the last dimension is index of network
varEstimComb = squeeze(cell2mat(varEstimComb));


%% A summary figure
figure
diffPosi = NetPars.Posi(3,:) - NetPars.Posi(1,:);
sumOAvgArray = squeeze(mean(OAvgArray, 1))';

diffMean = bsxfun(@minus, meanEstimNetSim(:, 2:end), ...
    meanEstimNetSim(:, 1));

subplot(2,2,1)
plot(diffPosi, meanEstimNetSim(:,1), 'linew', 1.5)
hold on
plot(diffPosi, meanEstimNetSim(:,2), 'r', 'linew', 1.5)
plot(diffPosi, meanEstimComb(:,1), 'k', 'linew', 1.5)
set(gca, 'xtick', -180: 90: 180, 'fontsize', 12)
xlabel('Cue disparity')
ylabel('Mean of S_1')
legend('Congruent', 'Opposite', 'Linear readout')
axis tight; axis square

subplot(2,2,2)
plot(diffPosi, diffMean, 'linew', 1.5)
hold on
% plot(diffPosi, diff(meanEstimComb, 1, 2), 'c', 'linew', 1.5)
plot([-180,180], [-180, 180], '--k')
axis(180*[-1 1 -1 1])
set(gca, 'xtick', -180: 90: 180, 'ytick', -180: 90: 180, ...
    'fontsize', 12)
legend('1O-1C', '2C-1C', '2O-1C', ...
    'location', 'northwest')
xlabel('Cue disparity'); ylabel('Estim disparity')
axis square

% subplot(2,2,2)
% plot(diffPosi, meanEstimComb(:,1), 'linew', 1.5)
% set(gca, 'xtick', -180: 90: 180, 'fontsize', 12)
% xlabel('Cue disparity'); ylabel('Mean of S_1')
% axis tight; axis square

subplot(2,2,3)
hold on
plot(diffPosi, sumOAvgArray(:, 1), 'linew', 1.5)
plot(diffPosi, sumOAvgArray(:, 2), 'r', 'linew', 1.5)
plot(diffPosi, sumOAvgArray(:, 3), '--', 'linew', 1.5)
plot(diffPosi, sumOAvgArray(:, 4), '--r', 'linew', 1.5)
axis tight; axis square
set(gca, 'xtick', -180: 90: 180, 'fontsize', 12)
xlabel('Cue disparity')
ylabel('Mean firing rate ')
legend('Net 1 congruent', 'Net 1 opposite', ...
    'Net 2 congruent', 'Net 2 opposite')

subplot(2,2,4)
hold on
weightComb = sumOAvgArray(:,1) ./ sum(sumOAvgArray(:,1:2), 2);

plot(diffPosi, weightComb.*varEstimNetSim(:, 1), 'linew', 1.5)
plot(diffPosi, (1-weightComb).*varEstimNetSim(:, 2), 'r', 'linew', 1.5)
plot(diffPosi, weightComb.*varEstimNetSim(:, 1) ...
    + (1-weightComb).*varEstimNetSim(:, 2), ...
    'k--', 'linew', 1.5)
plot(diffPosi, varEstimComb(:, 1), 'k', 'linew', 1.5)
axis tight; axis square
set(gca, 'xtick', -180: 90: 180, 'fontsize', 12)
xlabel('Cue disparity')
ylabel('Variance of  S_1')
legend('Congruent', 'Opposite', '')
% legend('Net 1 congruent', 'Net 1 opposite', ...
%     'Net 2 congruent', 'Net 2 opposite')


%%
figure
plotTuneCurveDemo