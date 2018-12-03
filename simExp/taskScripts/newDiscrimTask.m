% Estimate the neurometric function of a neuron (ROC analysis)

% Wen-Hao Zhang, Nov-1, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Parameters of inputs
% NetPars.AmplRatio = [0.4, 0.1, 0.5; 0.3, 0.9, 1.2];
NetPars.AmplRatio = [0.25, 0, 0.25; ...
                      0, 0.8, 0.8];
% NetPars.AmplRatio = [0.5, 0, 0.5; 0, 1.2, 1.2];
NetPars.AmplRatio = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet, 1]);
% NetPars.cueCond = [1,2,0]; % cue 1, cue 2, both cues
NetPars.cueCond = 0;

NetPars.tLen = 100 * NetPars.tau;
NetPars.nTrials_ROC = 30;

% Stimulus position
switch typeDiscrimTask
    case 1 
        % Discriminate whether x_1 is larger than 0 deg
        NetPars.Posi = repmat(-32:4:32, 2, 1);
        NetPars.Posi = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
        dAngle = NetPars.Posi(1,:);
    case 2
        % Discriminate whether x_1 is larger than x_2
        NetPars.Posi    = repmat(-16:2:16, 2, 1);
        NetPars.Posi(2,:) = flip(NetPars.Posi(2,:), 2);
        NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
        dAngle = NetPars.Posi(1,:) - NetPars.Posi(2,:);
end
   
% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);
parGrid = permute(parGrid, [2, 1]);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Net Simulation
RateNet = cell(size(parGrid));

parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/%d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    % Network input
    InputSet = makeNetInput([], netpars, struct('Iext', [], 'initRandStream', []));
    
    % Run simulation
    InputSet = simDecenNet(InputSet, netpars);
    RateNet{iterPar} = InputSet.O(:,:,netpars.tStat/netpars.dt+1:end, :);
end

RateNet = cellfun(@(x) reshape(x, 180, 4, [], NetPars.nTrials_ROC), RateNet, 'uniformout', 0);

%% ROC analysis
% find the most sensitive neuron
[~, IdxNeuron] = min(abs(NetPars.PrefStim - 90));

rateNeuron = cellfun(@(x) squeeze(x(IdxNeuron, [1,3],:,:)), RateNet, 'uniformout', 0);
rateNeuron = shiftdim(rateNeuron, -3);
rateNeuron = cell2mat(rateNeuron); % 5D array, [IdxGroup, Time, Trial, Posi, cueCond];
rateNeuron = permute(rateNeuron, [2,4,3,5,1]); % 5D array, [Time, Posi, Trial, cueCond, IdxGroup];

rateArray = cat(1, rateNeuron(:, (end+1)/2:end,:,:,:),...
    rateNeuron(:, (end+1)/2:-1:1,:,:,:)); % 2D matrix [2*Time, Posi * Trial * cueCond * numGroup]
sizeRate = size(rateArray);
rateArray = reshape(rateArray, sizeRate(1), []);

labels = [zeros(sizeRate(1)/2,1); ones(sizeRate(1)/2, 1)];

AUC = zeros(1, size(rateArray, 2)); % Area under ROC curve
for iter = 1: size(rateArray, 2)
    [~, ~, ~, AUC(iter)] = perfcurve(labels, rateArray(:, iter), 1);
end
AUC = reshape(AUC, sizeRate(2), []); % 2D array, [N, Trial*numPars*numGroup];
AUC = cat(1, 1 - flip(AUC(2:end,:), 1), AUC);

% Adjust the direction of AUC to make it into an overal ascending tendency
dAUC = sum(diff(AUC, 1, 1), 1);
Idx = (dAUC <0);
AUC(:, Idx) = flip(AUC(:, Idx), 1);

% fit with cumulative Gaussian function
ROCPars = zeros(2, size(AUC, 2));
for iter = 1: size(AUC, 2)
    ROCPars(:, iter) = lsNeuroMetric(dAngle, AUC(:, iter)');
end

AUC = reshape(AUC, [size(AUC,1), sizeRate(3:end)]); % 4D array, [nStim+1, Trial, cueCond, numGroup]
ROCPars = reshape(ROCPars, [2, sizeRate(3:end)]); % 4D array, [2, Trial, cueCond, numGroup]

%% Bayesian prediction
% random permute the order of trials
% ROCPars = ROCPars(:, randperm(NetPars.nTrials), :, :);
ROCPars(2,:,4,:) = sqrt(1./sum(1./ROCPars(2,:,1:2,:).^2, 3)); % [2, Trial, cueCond, numGroup]

% Statistical test of threshold in combined cue and prediction

switch typeDiscrimTask
    case 1
        % congruent cell
        ThetaNet = squeeze(ROCPars(2, :, 3, 1));
        ThetaBayes = squeeze(ROCPars(2, :, 4, 1));
    case 2
        % opposite cell
        ThetaNet = squeeze(ROCPars(2, :, 3, 2));
        ThetaBayes = squeeze(ROCPars(2, :, 4, 2));
end
% [p, h] = signrank(Theta, Theta1);
[h, p] = ttest2(ThetaNet, ThetaBayes);
clear Theta Theta1

%% Plot results
IdxTrial = randperm(NetPars.nTrials_ROC,1);
% IdxTrial = 20;
figure
for iter = 1: 6
    hAxe(iter) = subplot(2,3,iter);
    axis square;
    hold on
end
clear iter

cSpec = [44, 40, 187;
    255, 138, 55;
    42, 194, 64;
    42, 194, 64]/255;
sSpec = 'os^';

NetPars.cueCond = [1, 2, 0];

for iterGroup = 1: NetPars.numGroupPerNet
    IdxGroup = 1 + (iterGroup-1) * NetPars.numNets;
    
    % Plot a segment of tuning curve around zero
    axes(hAxe(1 + 3*(iterGroup-1)) )
    nStimPlot = 8;
    
    if exist('RateNet', 'var')
        % RateNet [N, numGroup, Time, Trial, cueCond];
        RatePlot = cellfun(@(x) squeeze(x(IdxNeuron, IdxGroup, :,:)), RateNet, 'uniformout', 0);
        RatePlot = shiftdim(RatePlot, -2);
        RatePlot = cell2mat(RatePlot); % 4D array, [Time, Trial, Posi, cueCond]
        RatePlot = squeeze(mean(RatePlot, 1)); % average over times
        
        meanFireRate = squeeze(mean(RatePlot, 1));
        stdFireRate = squeeze(std(RatePlot, 0, 1));
        clear RatePlot
                
        for iterPar = 1: length(NetPars.cueCond)
            switch typeDiscrimTask
                case 1
                    xValue = dAngle;
                case 2
                    %                     if NetPars.cueCond(iterPar) == 2
                    %                         xValue = NetPars.Posi(2,:);
                    %                     else
                    %                         xValue = NetPars.Posi(1,:);
                    %                     end
                    xValue = dAngle;
            end
            
            plot(xValue, ...
                meanFireRate(:, iterPar), ...
                [sSpec(iterPar), '-'], 'color', cSpec(iterPar,:), 'linew', 1);
            errorbar(xValue, ...
                meanFireRate(:, iterPar), ...
                stdFireRate(:, iterPar), ...
                [sSpec(iterPar), '-'], 'color', cSpec(iterPar,:), 'linestyle', 'none');
        end
    end
    
    % Plot an example of neurometric function
    axes(hAxe(2 + 3*(iterGroup-1)) )
    for iterPar = 1: length(NetPars.cueCond)
        plot(dAngle, mean(AUC(:, :, iterPar, iterGroup), 2), ...
            sSpec(iterPar), 'color', cSpec(iterPar,:), 'linew',1);
%         plot(dAngle, mean(AUC(:, IdxTrial, iterPar, iterGroup), 2), ...
%             sSpec(iterPar), 'color', cSpec(iterPar,:), 'linew',1);
        xAngle = linspace(dAngle(round(end/2)-nStimPlot), dAngle(round(end/2)+nStimPlot), 100);
        cdfPars = lsNeuroMetric(dAngle, mean(AUC(:,:,iterPar,iterGroup), 2)');
%         cdfPars = ROCPars(:, IdxTrial, iterPar, iterGroup);
        plot(xAngle, normcdf(xAngle, cdfPars(1), cdfPars(2)), ...
            'color', cSpec(iterPar,:), 'linew', 1);
    end
    
    % Plot the statistics of neuronal discrimination threshold
    axes(hAxe(3 + 3*(iterGroup-1)) )
    for iterPar = 1: 4
        bar(iterPar, mean(ROCPars(2, :, iterPar, iterGroup)), ...
            'facecolor', cSpec(iterPar,:))
    end
    errorbar(1:4, squeeze(mean(ROCPars(2,:,:,iterGroup),2)), ...
        squeeze(std(ROCPars(2,:,:,iterGroup),[],2)), ...
        'k', 'linew',1, 'linestyle', 'none')
    
    % Figure format
    axes(hAxe(1 + 3*(iterGroup-1)) )
    xlim([xValue(1)-2, xValue(end)+2])
    xlabel('Stimulus direction', 'fontsize', 9)
    ylabel('Firing rate (Hz)', 'fontsize', 9)
    set(gca, 'xtick', sort([dAngle([1, end]),0]) )
    
    axes(hAxe(2 + 3*(iterGroup-1)) )
    xlim([dAngle(1)-2, dAngle(end)+2])
    title(['Trial: ' num2str(IdxTrial)])
    switch typeDiscrimTask
    case 1
        xlabel('Stimulus direction', 'fontsize', 9)
    case 2
        xlabel('Disparity between two cues (x_1-x_2)', 'fontsize', 9)
    end
    ylabel('Correct fraction (%)', 'fontsize', 9)
    set(gca, 'xtick', sort([dAngle([1, end]),0]) )
    
    axes(hAxe(3 + 3*(iterGroup-1)) )
    xlim([0.5, 4.5])
    set(gca, 'xtick', 1:4, 'xticklabel', {'Cue 1','Cue 2','Combined cue','Prediction'}, ...
        'xticklabelrotation', 45)
    ylabel('Threshold', 'fontsize', 9)
    % ylim([10, 25])
end

switch typeDiscrimTask
    case 1
        axes(hAxe(3))
        title(['p=' num2str(p)])
    case 2
        axes(hAxe(6))
        title(['p=' num2str(p)])
end

% Specify input intensity
axes(hAxe(1))
title(['\alpha_1 = ' num2str(NetPars.AmplRatio(1)), ' , \alpha_2 =' num2str(NetPars.AmplRatio(2))]);


for iter = 1: length(hAxe)
    set(hAxe(iter), 'linew',1, 'fontsize', 9)
end