% Demonstration of network population activities and estimation restuls

% Wen-Hao Zhang, Apr-23, 2016

%% Net simulation
setPath;
% Load parameters
paramsNetDemoVM;

NetPars.tLen = 20 * NetPars.tau;
NetPars.Jrc = 0.4 * NetPars.Jc;
NetPars.JrpRatio = 0.9;
NetPars.Ampl = 0.2 * ones(2,1) * NetPars.Uc;

NetPars.NoiseFree = 0; % 1: no noise; 0: add noise;
NetPars.Posi = [0, 100]';

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

%% Produce random seeds
[seedIntNoisArray, seedExtNoisArray] = initRandSeed(NetPars, dimPar, size(parGrid));

%% Net Simulation
meanNetSim = cell(size(parGrid));
varNetSim = cell(size(parGrid));
concNetSim = cell(size(parGrid)); % concentration parameter of circular statistics

RateAvg = cell(size(parGrid));
RateT = cell(size(parGrid));

% parpool(6);
tStart = clock;
for iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedIntNois = seedIntNoisArray(iterPar);
    netpars.seedExtNois = seedExtNoisArray(iterPar);
    
    if netpars.Jrc*(1 + netpars.JrpRatio) > netpars.Jc
        % Parameters of persistent activity
        BumpPos = nan(1,2);
        meanBumpPos = nan(4,1);
        varBumpPos = nan(4,1);
        concBumpPos = nan(4,1);
    else
        % Network input
        InputSet = makeNetInputVM(netpars.Posi, netpars);
        
        % Run simulation
        InputSet = simDecenNet2VM(InputSet, netpars);
        
        % Bump position
%         O = InputSet.O(:,:, 50*NetPars.tau/ NetPars.dt+1 : end);
        O = InputSet.O(:,:, 2*NetPars.tau/ NetPars.dt+1: end);
        [BumpPos, meanBumpPos, varBumpPos, concBumpPos] ...
            = statBumpPos(O, NetPars);
        
        RateAvg{iterPar} = mean(O,3);
        RateT{iterPar} = squeeze(mean(InputSet.O,1));
    end
    meanNetSim{iterPar} = shiftdim(meanBumpPos, -ndims(parGrid));
    varNetSim{iterPar} = shiftdim(varBumpPos, -ndims(parGrid));
    concNetSim{iterPar} = shiftdim(concBumpPos, -ndims(parGrid));
end

meanNetSim = cell2mat(meanNetSim); % the last dimension is index of network
varNetSim = cell2mat(varNetSim);
concNetSim = cell2mat(concNetSim);

tEnd = clock;

%% Bayes prediction
bayesRes = getBayesPrediction(meanNetSim, varNetSim, dimPar);
meanNetSimRad = meanNetSim * pi/ NetPars.Width;
bayesResVM = getBayesPredictionVM(meanNetSimRad, concNetSim, NetPars, dimPar);
clear meanNetSimRad

%% Net estimation results
figure; 
sSpec = 'os+x';
cSpec = 'bbrr';
for iter = 1: 2
    hAxe(iter) = subplot(2,2,iter);
    hold on
end

for iterGroup = 1:4
    plot(hAxe(1), bayesResVM.meanNetBayes(iterGroup), meanNetSim(1,:,iterGroup), ...
        sSpec(iterGroup), 'color', cSpec(iterGroup), 'markersize', 8);
    plot(hAxe(2), bayesResVM.concNetBayes(iterGroup), concNetSim(1,:,iterGroup), ...
        sSpec(iterGroup), 'color', cSpec(iterGroup), 'markersize', 8);
end
plot(hAxe(1), [0, 180], [0, 180], '--k')
plot(hAxe(2), [0, 2e3], [0, 2e3], '--k')
axes(hAxe(1))
xlabel('Bayesian mean')
ylabel('Network mean')
axes(hAxe(2))
xlabel('Bayesian concentration')
ylabel('Network concentration')

hAxe(3) = subplot(2,2, [3,4]);
hold on
t = (1: length(RateT{1})) * NetPars.dt;
plot(hAxe(3), t, RateT{1}(1,:));
plot(hAxe(3), t, RateT{1}(3,:), 'r');
xlabel('Time (\tau)')
ylabel('Firing Rate')
clear sumO

%% Figure
% 2D firing rate
IdxGroup = [1,3]; % 1: congruent neurons in net 1
% 3: opposite neurons in net 1
cSpec = [0,0,1;0,0,1; 1,0,0; 1,0,0];
titleCond = {'Combined','Cue 1','Cue 2'};

rPlot = 10:10:40;
thetaPlot = 0: 45: 135;

figure;
for iter = 1: 3
    hAxe(iter) = subplot(1,3,iter);
    hold on; axis square
end
for iterCond = 1:3
    axes(hAxe(iterCond));
    title(titleCond{iterCond});
    % Guide lines
    for iter = 1: length(rPlot)-1
        hPlot = polar(linspace(0, 2*pi, 360), rPlot(iter)* ones(1, 360), '--k');
    end
    hPlot = polar(linspace(0, 2*pi, 360), rPlot(end)* ones(1, 360), 'k');
    
    for iter = 1: length(thetaPlot)
        polar(thetaPlot(iter)*pi/180*ones(1,2), rPlot(end)*[-1, 1], '--k');
    end
    
    for iterGroup = IdxGroup
        ratePlot = [RateAvg{iterCond}(:, iterGroup); RateAvg{iterCond}(1, iterGroup)];
        
        hPlot = polar([NetPars.PrefStim; NetPars.PrefStim(1)]*pi/180, ratePlot);
        set(hPlot, 'color', cSpec(iterGroup,:))
        
        % Bump position
        polar(meanNetSim(iterCond,1, iterGroup)*pi/ 180 * ones(1,2), ...
            [0, max(RateAvg{iterCond}(:, iterGroup))], '--k')
    end
    axis square;axis tight;axis off;
    axis([-1 1 -1 1]*rPlot(end))
end

axes(hAxe(1))
% Cue position
polar(NetPars.Posi(1)* pi/ 180 * ones(1,2), ...
    [0, max(RateAvg{2}(:, iterGroup))], '--k')
polar(NetPars.Posi(2)* pi/ 180 * ones(1,2), ...
    [0, max(RateAvg{2}(:, iterGroup))], '--k')

%% Figure
% 3D firing rate
IdxGroup = [1,3]; % 1: congruent neurons in net 1
% 3: opposite neurons in net 1
cSpec = [0,0,1;0,0,1; 1,0,0; 1,0,0];
titleCond = {'Combined','Cue 1','Cue 2'};

figure;
for iter = 1: 3
    hAxe(iter) = subplot(1,3,iter);
    hold on; axis square
end

xPlot = cos(NetPars.PrefStim* pi/ 180);
yPlot = sin(NetPars.PrefStim* pi/ 180);
xPlot = [xPlot; xPlot(1)];
yPlot = [yPlot; yPlot(1)];

for iterCond = 1:3
    axes(hAxe(iterCond));
    title(titleCond{iterCond});
    for iterGroup = IdxGroup
        plot(0,0, 'ok')
        plot([0, cos(NetPars.Posi(1)* pi/ 180)], ...
            [0, sin(NetPars.Posi(1)* pi/ 180)], ...
            '--k');
        plot([0, cos(NetPars.Posi(2)* pi/ 180)], ...
            [0, sin(NetPars.Posi(2)* pi/ 180)], ...
            '--k');
        plot([0, -cos(NetPars.Posi(2)* pi/ 180)], ...
            [0, -sin(NetPars.Posi(2)* pi/ 180)], ...
            '--k');
        plot(xPlot,yPlot, 'k');
        
        plot3(xPlot,yPlot,  [RateAvg{iterCond}(:, iterGroup); RateAvg{iterCond}(1, iterGroup)], ...
            'color', cSpec(iterGroup,:));
        
        % Bump position
        xBumpPos = cos(meanNetSim(iterCond,1, iterGroup)*pi/ 180);
        yBumpPos = sin(meanNetSim(iterCond,1, iterGroup)*pi/ 180);
        plot3(xBumpPos*ones(1,2), yBumpPos*ones(1,2), [0, max(RateAvg{iterCond}(:, iterGroup))], 'k')
        
        view([-0, 20])
        axis off
    end
    zlim([0, 20])
end