% Demo of the flexible integration/segregation network model
% Wen-Hao Zhang, May-7, 2017
% wenhaoz1@andrew.cmu.edu
% Carnegie Mellon University

setWorkPath;

% Parameters of inputs
NetPars.Posi = repmat([-30; 30], [NetPars.numGroupPerNet, 1]);
% NetPars.Posi = repmat([0; 0], [NetPars.numGroupPerNet, 1]);

NetPars.AmplRatio = 0.35; %0.7
NetPars.tLen      = 550 * NetPars.tau;
NetPars.cueCond   = [1,2,0];

% NetPars.tLen    = 20 * NetPars.tau;
% NetPars.nTrials = 50;
% NetPars.tStat   = 15 * NetPars.tau; % The starting time to make statistics

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);
%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Simulation of Decision-making circuit
NetStatStruct = struct('BumpPos', [], ...
    'OHeight', [], ...
    'meanBumpPos', [], ...
    'mrlBumpPos', [], ...
    'concBumpPos', [], ...
    'varBumpPos', [], ...
    'OHeightAvg', [], ...
    'OAvgXTime', [], ...
    'OStdXTime', []);
NetStatStruct = repmat(NetStatStruct, size(parGrid));
NetRes = cell(size(parGrid));

tStart = clock;
for iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/%d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    % Network input
    InputSet = makeNetInput([], netpars, struct('Iext', [], 'initRandStream', []));
    %     InputSet.Iext = repmat(InputSet.Iext, [1,1, InputSet.szIext(3)]);
    %     InputSet.Iext(:,:,1:2e3) = 0;
    %     InputSet.Iext(:,:,1.2e4+1:end) = 0;
    
    % Run simulation
    outArgs = struct('InputSet', [], ...
        'NetStat', NetStatStruct(iterPar));
    [InputSet, NetStatStruct(iterPar)] = simDecenNet(InputSet, netpars, outArgs);
    NetRes{iterPar} = InputSet.O;
end
tEnd = clock;
clear netpars

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

%% Plot population activities under three cueing conditions
IdxGroup = 1;
IdxCueCond = 3;
cSpec = lines(2);

figure
subplot(3,3,1:2)
tMax = 1.4e4;
tPlot = (1:tMax) * NetPars.dt/NetPars.tau;
colormap(jet)
imagesc(tPlot, NetPars.PrefStim, squeeze(NetRes{IdxCueCond}(:,1,1:tMax)))
axis xy
set(gca, 'ytick', [-178, -90:90:180], 'yticklabel', -180:90:180, ...
    'xtick', 0:50:150)
ylabel('Neuron index \theta')
xlabel('Time (\tau)')

subplot(3,3,3)
hold on
fill([NetEstimRes.OAvgXTime(IdxCueCond, :, IdxGroup)' + NetEstimRes.OStdXTime(IdxCueCond, :, IdxGroup)'; ...
    flipud(NetEstimRes.OAvgXTime(IdxCueCond, :, IdxGroup)' - NetEstimRes.OStdXTime(IdxCueCond, :, IdxGroup)')], ...
    [NetPars.PrefStim; flipud(NetPars.PrefStim)], ...
    lines(1), 'linestyle', '-', 'facecolor', cSpec(IdxGroup,:), ...
    'facealpha', 0.3)
plot(NetEstimRes.OAvgXTime(3, :, 1), NetPars.PrefStim, ...
    'color', cSpec(1,:), 'linew', 1.5);
set(gca, 'ytick', [-178, -90:90:180], 'yticklabel', -180:90:180, 'ylim', NetPars.PrefStim([1, end]));
xlabel('Firing rate (Hz)')

for iter = 1: 6
    hAxe(iter) = subplot(3,3, iter+3);
end

% IdxGroup = [1,3];
IdxGroup = 1:4;

% maxYVal = 1.2*max(NetEstimRes.OAvgXTime(:));
maxYVal = 60;
for iterCueCond = 1: 3
    for iterGroup = 1: length(IdxGroup)
        axes(hAxe(iterCueCond+(1-mod(iterGroup, 2)) * 3)); hold on
        
        idxGroup = IdxGroup(iterGroup);
        plot(NetPars.PrefStim, NetEstimRes.OAvgXTime(iterCueCond, :, idxGroup), ...
            'color', cSpec(round(iterGroup/2),:), 'linew', 1.5);
        fill([NetPars.PrefStim; flipud(NetPars.PrefStim)], ...
            [NetEstimRes.OAvgXTime(iterCueCond, :, idxGroup)' + NetEstimRes.OStdXTime(iterCueCond, :, idxGroup)'; ...
            flipud(NetEstimRes.OAvgXTime(iterCueCond, :, idxGroup)' - NetEstimRes.OStdXTime(iterCueCond, :, idxGroup)')], ...
            lines(1), 'linestyle', '-', 'facecolor', cSpec(round(iterGroup/2),:), ...
            'facealpha', 0.3)
        
        % Plot the dashed line indicating bump position
        plot(NetEstimRes.meanBumpPos(iterCueCond, idxGroup) * ones(1,2), ...
            [0, maxYVal], '--', 'color', cSpec(round(iterGroup/2),:));
        switch iterCueCond
            case 1
                plot(NetPars.Posi(1) * ones(1,2), [0, maxYVal], '--k');
            case 2
                plot(NetPars.Posi(2) * ones(1,2), [0, maxYVal], '--k');
            case 3
                plot(NetPars.Posi(2 - mod(iterGroup,2)) * ones(1,2), [0, maxYVal], '--k');
        end
        
        set(gca, 'xlim', NetPars.PrefStim([1, end]), 'xtick', [-178, -90:90:180], ...
            'xticklabel', [-180, -90:90:180], 'ylim', [0, maxYVal], 'ytick', 0:10:maxYVal)
        %     axis square
    end
end
axes(hAxe(2))
xlabel('Neuron index (\theta)')
axes(hAxe(1))
ylabel('Firing rate (Hz)')
axes(hAxe(4))
title(['\alpha_1=\alpha_2=' num2str(NetPars.AmplRatio)])

%% Plot the network estimate (bump position) at three cueing conditions
% Joint posterior under three cueing conditions

tStep = 1e2;
alphaValue = 1;
cSpec = [44, 40, 187; % cue 1
    255, 138, 55; % cue 2
    42, 194, 64;]/255; % combined cues

IdxGroup = [1,2; 3,4];
% Iterate on congruent and opposite cells
for iterGroup = 1: size(IdxGroup, 1)
    % Generate new figure
    figure;
    hAxe(1) = axes('position', [0.1, 0.1, 0.6, 0.6]); hold on; daspect([1,1,1])
    hAxe(2) = axes('position', [0.1, 0.75, 0.6, 0.15]); hold on; %daspect([1,4,1])
    hAxe(3) = axes('Position', [0.75, 0.1, 0.15, 0.6]); hold on; %daspect([4,1,1])
    
    BumpPos = NetEstimRes.BumpPos(:,IdxGroup(iterGroup,:), 1e4:tStep:NetPars.tLen/NetPars.dt);
    
    % Determine the boundary of axis
    axisLim = [min(reshape(BumpPos(:,1,:),1,[])), ...
        max(reshape(BumpPos(:,1,:),1,[])), ...
        min(reshape(BumpPos(:,2,:),1,[])), ...
        max(reshape(BumpPos(:,2,:),1,[]))];
    axisLim(1:2) = axisLim(1:2) + 0.1*diff(axisLim(1:2)) * [-1,1];
    axisLim(3:4) = axisLim(3:4) + 0.1*diff(axisLim(3:4)) * [-1,1];
    axisLim = round(axisLim/10) * 10;
    switch iterGroup
        case 1  % congruent cell
            axisLim(1:2) = max(abs(axisLim(1:2))) * [-1, 1];
            axisLim(3:4) = max(abs(axisLim(3:4))) * [-1, 1];
        case 2  % opposite cell
            axisLim([1,4]) = max(abs(axisLim([1,4]))) * [-1, 1];
            axisLim([2,3]) = max(abs(axisLim([2,3]))) * [-1, 1];
    end

%     cSpec = lines(3);
    for iterCueCond = 1:3
        scatter(hAxe(1), BumpPos(iterCueCond,1,:), BumpPos(iterCueCond,2,:), ...
            'o', 'markerfacecolor', cSpec(iterCueCond,:), 'markeredgecolor', ones(1,3), ...
            'sizedata', 10);
        XAvg = squeeze(mean(BumpPos(iterCueCond,1:2,:),3));
        XCov = cov(squeeze(BumpPos(iterCueCond,1:2,:))');
       
        fh = @(x,y) ( ([x,y] - XAvg) / XCov/9* ([x,y]-XAvg)' - 1);
        hEllipse = ezplot(hAxe(1), fh, [XAvg(1) + 3*XCov(1)*[-1, 1], XAvg(2) + 3*XCov(4)*[-1, 1]]);
        
        cSpecA = (1-alphaValue)*ones(1,3) + alphaValue*cSpec(iterCueCond,:);
        set(hEllipse, 'color', cSpecA, 'linew', 1.5)
        plot(hAxe(1), XAvg(1),XAvg(2), 'x', 'linew', 2, 'color', cSpec(iterCueCond,:));
        
        % Marginal distribution of network 1's estimate
        histEdge = linspace(axisLim(1), axisLim(2), 2e2);
        hBumpPos = histc(reshape(BumpPos(iterCueCond,1,:),1,[]),histEdge);
        hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
        stairs(hAxe(2), histEdge, hBumpPos, 'color', cSpec(iterCueCond,:));
        plot(hAxe(2), histEdge, normpdf(histEdge, XAvg(1), sqrt(XCov(1))), 'color', cSpecA, ...
            'linew', 2);
        
        % Marginal distribution of network 2's estimate
        histEdge = linspace(axisLim(3), axisLim(4), 2e2);
        hBumpPos = histc(reshape(BumpPos(iterCueCond,2,:),1,[]),histEdge);
        hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
        stairs(hAxe(3), hBumpPos, histEdge, 'color', cSpec(iterCueCond,:));
        plot(hAxe(3), normpdf(histEdge, XAvg(2), sqrt(XCov(2,2))), histEdge, 'color', cSpecA, ...
            'linew', 2);
    end
    
    axes(hAxe(1))
    switch iterGroup
        case 1
            % Congruent neuron
            plot(NetPars.Posi(1)*ones(1,2), axisLim(3:4), '--k')
            plot(NetPars.Posi(2)*ones(1,2), axisLim(3:4), '--k')
            plot(axisLim(1:2), NetPars.Posi(2)*ones(1,2), '--k')
            plot(axisLim(1:2), NetPars.Posi(1)*ones(1,2), '--k')
            plot(NetPars.Posi(1:2), NetPars.Posi(1:2), '--k')
            plot(NetPars.Posi(1:2), NetPars.Posi([2,1]), '--k')
            
            plot(NetPars.Posi(1:2), NetPars.Posi(1:2), '--k')
            plot(NetPars.Posi(1:2), NetPars.Posi([2,1]), '--k')
            
            xlabel('Congruent neurons in module 1, p(s_1|x_1,x_2)')
            ylabel('Congruent neurons in module 2, p(s_2|x_1,x_2)')
            set(hAxe(1), 'xtick', NetPars.Posi(1:2), 'ytick', NetPars.Posi(1:2))
        case 2
            % Opposite neuron
            PosiOppo = angle(exp(1i*NetPars.Posi(1:2)*pi/NetPars.Width+pi*1i))*NetPars.Width/pi;
            plot(NetPars.Posi(1)*ones(1,2), axisLim(3:4), '--k')
            plot(PosiOppo(2)*ones(1,2), axisLim(3:4), '--k')
            plot(axisLim(1:2), NetPars.Posi(2)*ones(1,2), '--k')
            plot(axisLim(1:2), PosiOppo(1)*ones(1,2), '--k')
            
            plot([NetPars.Posi(1), PosiOppo(2)], [NetPars.Posi(2), PosiOppo(1)], '--k')
            plot([NetPars.Posi(1), PosiOppo(2)], [PosiOppo(1), NetPars.Posi(2)], '--k')
            
            xlabel('Opposite neurons in module 1, D(s_1|x_1,x_2)')
            ylabel('Opposite neurons in module 2, D(s_2|x_1,x_2)')
            set(hAxe(1), 'xtick', sort([NetPars.Posi(1), PosiOppo(2)], 'ascend'), ...
                'ytick', [NetPars.Posi(2), PosiOppo(1)])
    end
    axis(axisLim)
    axis square; axis tight
    box on
    title('')
    
    set(hAxe(2:3), 'xtick', {}, 'ytick', {}, 'xticklabel', {}, 'yticklabel', {})
    
    linkaxes(hAxe(1:2), 'x')
    linkaxes(hAxe([1,3]), 'y')
    axes(hAxe(2)); axis tight
    axes(hAxe(3)); axis tight
end

%%
plotJointDistribution_170619;

%% Demonstration of reconstruction the estimates of stimuli under single cues.

DemoDecode_SingleCues;

