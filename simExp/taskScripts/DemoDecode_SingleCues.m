% Demonstration of decoding cues from congruent and opposite neurons

% Wen-Hao Zhang, June-22, 2017
% wenhaoz1@andrew.cmu.edu
% Carnegie Mellon University

flagDecodeMethod = 2;

IdxCombCue = find(NetPars.cueCond==0);

switch flagDecodeMethod
    case 1
        vecNetRes = NetEstimRes.OHeight(IdxCombCue,:,:) .* ...
            exp(1i*NetEstimRes.BumpPos(IdxCombCue,:,:)*pi/NetPars.Width);

        % size of decodeNetRes: [cue, network, time]
        decodeNetRes = cat(2, ...
            (vecNetRes(:,1,:) + vecNetRes(:,3,:))/2, ... % p(s_1|x_1)
            (vecNetRes(:,1,:) - vecNetRes(:,3,:))/2, ... % p(s_1|x_2)
            (vecNetRes(:,2,:) - vecNetRes(:,4,:))/2, ... % p(s_2|x_1)
            (vecNetRes(:,2,:) + vecNetRes(:,4,:))/2);% p(s_2|x_2)
        decodeNetRes = reshape(decodeNetRes, 2, 2, []);
        
        meanDecodeNet = angle(decodeNetRes) * NetPars.Width/pi;
        concDecodeNet = abs(decodeNetRes);
    case 2
        ReconstPopRes = cat(2, ...
            sum(NetRes{IdxCombCue}(:,[1,3],:),2), ... % p(s_1|x_1)
            diff(NetRes{IdxCombCue}(:,[3,1],:),1, 2), ... % p(s_1|x_2)
            diff(NetRes{IdxCombCue}(:,[4,2],:),1, 2), ... % p(s_2|x_1)
            sum(NetRes{IdxCombCue}(:,[2,4],:),2)); % p(s_2|x_2)
        
        meanDecodeNet = statBumpPos(ReconstPopRes, NetPars);
        meanDecodeNet = reshape(meanDecodeNet, 2, 2, []);
end

%% Joint distribution
figure
hAxe(1) = axes('position', [0.1, 0.1, 0.6, 0.6]); hold on; daspect([1,1,1])
hAxe(2) = axes('position', [0.1, 0.75, 0.6, 0.15]); hold on; %daspect([1,4,1])
hAxe(3) = axes('Position', [0.75, 0.1, 0.15, 0.6]); hold on; %daspect([4,1,1])

tStep = 100;

cSpec = [44, 40, 187; % cue 1
    255, 138, 55; % cue 2
    42, 194, 64;]/255; % combined cues
IdxCueCond = 3;
alphaValue = 1;

% Congruent and opposite cells under combined cues
for iterGroup = 1: size(IdxGroup, 1)
    BumpPos = NetEstimRes.BumpPos(:,IdxGroup(iterGroup,:), 1e4:tStep:NetPars.tLen/NetPars.dt);
    
    if iterGroup == 2
        cSpec = 1 - cSpec;
    end
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
    
    
    plot(hAxe(1), squeeze(BumpPos(IdxCueCond,1,:)), squeeze(BumpPos(IdxCueCond,2,:)), ...
        '.', 'color', cSpec(IdxCueCond,:));
    %     scatter(hAxe(1), BumpPos(IdxCueCond,1,:), BumpPos(IdxCueCond,2,:), ...
    %         'o', 'markerfacecolor', cSpec(IdxCueCond,:), 'markeredgecolor', ones(1,3), ...
    %         'sizedata', 10);
    XAvg = squeeze(mean(BumpPos(IdxCueCond,1:2,:),3));
    XCov = cov(squeeze(BumpPos(IdxCueCond,1:2,:))');
    
    fh = @(x,y) ( ([x,y] - XAvg) / XCov/9* ([x,y]-XAvg)' - 1);
    hEllipse = ezplot(hAxe(1), fh, [XAvg(1) + 3*XCov(1)*[-1, 1], XAvg(2) + 3*XCov(4)*[-1, 1]]);
    
    cSpecA = (1-alphaValue)*ones(1,3) + alphaValue*cSpec(IdxCueCond,:);
    set(hEllipse, 'color', cSpecA, 'linew', 1.5)
    plot(hAxe(1), XAvg(1),XAvg(2), 'x', 'linew', 2, 'color', cSpec(IdxCueCond,:));
    
    % Marginal distribution of network 1's estimate
    histEdge = linspace(axisLim(1), axisLim(2), 2e2);
    hBumpPos = histc(reshape(BumpPos(IdxCueCond,1,:),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(2), histEdge, hBumpPos, 'color', cSpec(IdxCueCond,:));
    plot(hAxe(2), histEdge, normpdf(histEdge, XAvg(1), sqrt(XCov(1))), 'color', cSpecA, ...
        'linew', 2);
    
    % Marginal distribution of network 2's estimate
    histEdge = linspace(axisLim(3), axisLim(4), 2e2);
    hBumpPos = histc(reshape(BumpPos(IdxCueCond,2,:),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(3), hBumpPos, histEdge, 'color', cSpec(IdxCueCond,:));
    plot(hAxe(3), normpdf(histEdge, XAvg(2), sqrt(XCov(2,2))), histEdge, 'color', cSpecA, ...
        'linew', 2);
    
    axes(hAxe(1))
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

% ------------------------------------------------------------------------
% Decoding resutls
cSpec = [44, 40, 187; % cue 1
    255, 138, 55; % cue 2
    42, 194, 64;]/255; % combined cues

alphaValue = 0.5;

for iterCueCond = 1:2
    BumpPos = NetEstimRes.BumpPos(:,1:2,1e4+1:tStep:end);
    % Determine the boundary of axis
    axisLim = [min(reshape(BumpPos(:,1,:),1,[])), ...
        max(reshape(BumpPos(:,1,:),1,[])), ...
        min(reshape(BumpPos(:,2,:),1,[])), ...
        max(reshape(BumpPos(:,2,:),1,[]))];
    axisLim(1:2) = axisLim(1:2) + 0.1*diff(axisLim(1:2)) * [-1,1];
    axisLim(3:4) = axisLim(3:4) + 0.1*diff(axisLim(3:4)) * [-1,1];
    axisLim = round(axisLim/10) * 10;
    
    % p(s_1,s_2|x_iterCueCond)
    % ---------------------------------------------------------------------
    % Actual results
    cSpecActual = (1-alphaValue)*ones(1,3) + alphaValue*cSpec(iterCueCond,:);
    plot(hAxe(1), squeeze(BumpPos(iterCueCond,1,:)), ...
        squeeze(BumpPos(iterCueCond,2,:)), '.', 'color', ...
        (1-alphaValue)*ones(1,3) + alphaValue*cSpec(iterCueCond,:));
    
    XAvg = squeeze(mean(BumpPos(iterCueCond,1:2,:),3));
    XCov = cov(squeeze(BumpPos(iterCueCond,1:2,:))');
    fh = @(x,y) ( ([x,y] - XAvg) / XCov/9* ([x,y]-XAvg)' - 1);
    hEllipse = ezplot(hAxe(1), fh, [XAvg(1) + 3*XCov(1)*[-1, 1], XAvg(2) + 3*XCov(4)*[-1, 1]]);
    set(hEllipse, 'color', cSpecActual, 'linew', 1.5)
    
    % Marginal distribution of network 1's estimate
    histEdge = linspace(axisLim(1), axisLim(2), 2e2);
    hBumpPos = histc(reshape(BumpPos(iterCueCond,1,:),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(2), histEdge, hBumpPos, 'color', cSpecActual);
    plot(hAxe(2), histEdge, normpdf(histEdge, XAvg(1), sqrt(XCov(1))), 'color', cSpecActual, ...
        'linew', 2);
    
    % Marginal distribution of network 2's estimate
    histEdge = linspace(axisLim(3), axisLim(4), 2e2);
    hBumpPos = histc(reshape(BumpPos(iterCueCond,2,:),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(3), hBumpPos, histEdge, 'color', cSpecActual);
    plot(hAxe(3), normpdf(histEdge, XAvg(2), sqrt(XCov(2,2))), histEdge, 'color', cSpecActual, ...
        'linew', 2);
    
    % ---------------------------------------------------------------------
    % Reconstruction    
    plot(hAxe(1), squeeze(meanDecodeNet(iterCueCond,1,1e4+1:tStep:end)), ...
        squeeze(meanDecodeNet(iterCueCond,2,1e4+1:tStep:end)), '.', 'color', ...
        cSpec(iterCueCond,:))
    XAvg = squeeze(mean(meanDecodeNet(iterCueCond,1:2,1e4+1:tStep:end),3));
    XCov = cov(squeeze(meanDecodeNet(iterCueCond,1:2,1e4+1:tStep:end))');
    fh = @(x,y) ( ([x,y] - XAvg) / XCov/9* ([x,y]-XAvg)' - 1);
    hEllipse = ezplot(hAxe(1), fh, [XAvg(1) + 3*XCov(1)*[-1, 1], XAvg(2) + 3*XCov(4)*[-1, 1]]);
    set(hEllipse, 'color', cSpec(iterCueCond,:), 'linew', 1.5)
    
    
    % Marginal distribution of Reconstruction of network 1's estimate
    histEdge = linspace(axisLim(1), axisLim(2), 2e2);
    hBumpPos = histc(reshape(meanDecodeNet(iterCueCond,1,1e4+1:tStep:end),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(2), histEdge, hBumpPos, 'color', cSpec(iterCueCond,:));
    plot(hAxe(2), histEdge, normpdf(histEdge, XAvg(1), sqrt(XCov(1))), 'color', cSpec(iterCueCond,:), ...
        'linew', 2);
    
    % Marginal distribution of Reconstruction of network 2's estimate
    histEdge = linspace(axisLim(3), axisLim(4), 2e2);
    hBumpPos = histc(reshape(meanDecodeNet(iterCueCond,2,1e4+1:tStep:end),1,[]),histEdge);
    hBumpPos = hBumpPos/ (sum(hBumpPos)*mean(diff(histEdge)));
    stairs(hAxe(3), hBumpPos, histEdge, 'color', cSpec(iterCueCond,:));
    plot(hAxe(3), normpdf(histEdge, XAvg(2), sqrt(XCov(2,2))), histEdge, 'color', cSpec(iterCueCond,:), ...
        'linew', 2);
    
end


% ---------- Plot dashed lines -----------
% Determine the boundary of axis
BumpPos = NetEstimRes.BumpPos(1:2,1:2, 1e4:tStep:NetPars.tLen/NetPars.dt);
axisLim = [min(reshape(BumpPos(:,1,:),1,[])), ...
    max(reshape(BumpPos(:,1,:),1,[])), ...
    min(reshape(BumpPos(:,2,:),1,[])), ...
    max(reshape(BumpPos(:,2,:),1,[]))];
axisLim(1:2) = axisLim(1:2) + 0.1*diff(axisLim(1:2)) * [-1,1];
axisLim(3:4) = axisLim(3:4) + 0.1*diff(axisLim(3:4)) * [-1,1];
axisLim = round(axisLim/10) * 10;

plot(hAxe(1), axisLim(1:2), NetPars.Posi(1,:)*ones(1,2), '--k')
plot(hAxe(1), axisLim(1:2), NetPars.Posi(2,:)*ones(1,2), '--k')
plot(hAxe(1), NetPars.Posi(1,:)*ones(1,2), axisLim(3:4), '--k')
plot(hAxe(1), NetPars.Posi(2,:)*ones(1,2), axisLim(3:4), '--k')

axis(axisLim)
linkaxes(hAxe(1:2), 'x')
linkaxes(hAxe([1,3]), 'y')
axes(hAxe(2)); axis tight
axes(hAxe(3)); axis tight

axis([-60, 60, -60, 60]);
axes(hAxe(1))
title('')
xlabel('p(s_1|x_1,x_2)')
ylabel('p(s_2|x_1,x_2)')

%% Comparison between decoding and actual
figure
cSpec = lines(2);
hold on
% Direct cue of two networks
plot(squeeze(NetEstimRes.BumpPos(1,1,1e4+1:tStep:end)), ...
    squeeze(meanDecodeNet(1,1,1e4+1:tStep:end)), '.', 'color', cSpec(1,:))
% plot(squeeze(NetEstimRes.BumpPos(2,2,1e4+1:tStep:end)), ...
%     squeeze(meanDecodeNet(2,2,1e4+1:tStep:end)), '.', 'color', cSpec(1,:))

% Indirect cue of two networks
plot(squeeze(NetEstimRes.BumpPos(2,1,1e4+1:tStep:end)), ...
    squeeze(meanDecodeNet(2,1,1e4+1:tStep:end)), '.', 'color', cSpec(2,:))
% plot(squeeze(NetEstimRes.BumpPos(1,2,1e4+1:tStep:end)), ...
%     squeeze(meanDecodeNet(1,2,1e4+1:tStep:end)), '.', 'color', cSpec(2,:))

plot([-40, 60], [-40, 60], '--k')
axis square
xlabel('Actual mean of estimate of s_1')
ylabel('Reconstructed mean of estimate of s_1')
set(gca, 'xtick', -40:20:60, 'ytick', -40:20:60)
