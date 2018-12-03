% Generate new figure
figure;

tStep = 1e2;
alphaValue = 1;
cSpec = [44, 40, 187; % cue 1
    255, 138, 55; % cue 2
    42, 194, 64;]/255; % combined cues

IdxGroup = [1,2; 3,4];
hAxe(1) = axes('position', [0.1, 0.1, 0.6, 0.6]); hold on; daspect([1,1,1])
hAxe(2) = axes('position', [0.1, 0.75, 0.6, 0.15]); hold on; %daspect([1,4,1])
hAxe(3) = axes('Position', [0.75, 0.1, 0.15, 0.6]); hold on; %daspect([4,1,1])

% Iterate on congruent and opposite cells
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