%%
% Parameters for plots
% tStat = NetPars.tStat;
tStat = 6e3;

%%
axisBound = zeros(4,2,3);
for IdxCueCond = 1: length(NetPars.cueCond)
    for IdxGroup = 1: NetPars.numNets*NetPars.numGroupPerNet
        axisBound(IdxGroup,1,IdxCueCond) = min(reshape(NetEstimRes.BumpPos(IdxCueCond, :, tStat:end, IdxGroup), [], 1));
        axisBound(IdxGroup,2,IdxCueCond) = max(reshape(NetEstimRes.BumpPos(IdxCueCond, :, tStat:end, IdxGroup), [], 1));
    end
end

% tPlot = (1: NetPars.tTrial/NetPars.dt) * NetPars.dt/ NetPars.tau;
tPlot = (1: NetPars.tLen/NetPars.dt) * NetPars.dt/ NetPars.tau;
for iterTrial = 1: NetPars.nTrials
    figure(1)
    for IdxCueCond = 1: length(NetPars.cueCond)
        subplot(3,4, (IdxCueCond-1)*4+1);
        plot(tPlot, squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, :, 1)))
        hold on
        plot(tPlot, squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, :, 3)))
        axis square; hold off
        %     ylim([-40, 40])
        title(sprintf('Trial num.: %d.', iterTrial))
        xlim([tPlot(1), tPlot(end)])
        ylabel('Bump position (Net 1)')
        xlabel('Time (\tau)')
        
        subplot(3,4,(IdxCueCond-1)*4+2);
        plot(tPlot, squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, :, 1)))
        hold on
        plot(tPlot, squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, :, 3)))
        axis square; hold off
        xlim([tPlot(1), tPlot(end)])
        ylabel('Bump height')
        xlabel('Time (\tau)')
        
        subplot(3,4,(IdxCueCond-1)*4+3);
        plot(tPlot, squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, :, 2)))
        hold on
        plot(tPlot, squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, :, 4)))
        axis square; hold off
        %     ylim(90+ [-40, 40])
        xlim([tPlot(1), tPlot(end)])
        ylabel('Bump position (Net 2)')
        xlabel('Time (\tau)')
        
        subplot(3,4,(IdxCueCond-1)*4+4);
        plot(tPlot, squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, :, 2)))
        hold on
        plot(tPlot, squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, :, 4)))
        axis square; hold off
        xlim([tPlot(1), tPlot(end)])
        ylabel('Bump height')
        xlabel('Time (\tau)')
    end
    
    figure(2)
    strTitle = {'Combined', 'Cue 1', 'Cue 2'};
    
    for IdxCueCond = 1: length(NetPars.cueCond)
        subplot(3,2, 2*IdxCueCond-1)
        plot(squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:end, 1)), ...
            squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:end, 2)), '.')
        axis square
        xlabel('C neurons in net1')
        ylabel('C neurons in net2')
        title(strTitle{IdxCueCond})
        % axis([axisBound(1,1:2,IdxCueCond), axisBound(2,1:2,IdxCueCond)]);
        axis(NetPars.Posi')

        subplot(3,2,2*IdxCueCond)
        plot(squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:end, 3)), ...
            squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:end, 4)), '.')        
        axis square
        xlabel('O neurons in net1')
        ylabel('O neurons in net2')
        title(strTitle{IdxCueCond})
        axis([axisBound(3,1:2,IdxCueCond), axisBound(4,1:2,IdxCueCond)]);
    end
    
    pause
end


%% Scatter plot according to the winner 
figure
hold on
for iter = 1: 2
    hAxe(iter) = subplot(1,2,iter);
    hold on
end
clear iter

colorPan = lines(2);
IdxCueCond = 1;
for iterTrial = 1: NetPars.nTrials
    
    flagCWin = (mean(squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, tStat:end, 1))) > ...
        mean(squeeze(NetEstimRes.OHeight(IdxCueCond, iterTrial, tStat:end, 3))));
    if flagCWin
        cSpec = colorPan(1,:);
    else
        cSpec = colorPan(2,:);
    end
    
    plot(hAxe(1), squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:10:end, 1)), ...
        squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:10:end, 2)), '.', 'color', cSpec(1,:))
    plot(hAxe(2), squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:10:end, 3)), ...
        squeeze(NetEstimRes.BumpPos(IdxCueCond, iterTrial, tStat:10:end, 4)), '.', 'color', cSpec(1,:))
end
axes(hAxe(1))
axis square
axis(NetPars.Posi')
plot(NetPars.Posi(1:2)', NetPars.Posi(1:2)', 'k--')
plot(NetPars.Posi(1:2)', NetPars.Posi([2,1])', 'k--')
xlabel('C neurons in net1')
ylabel('C neurons in net2')
axes(hAxe(2))
axis square
xlabel('O neurons in net1')
ylabel('O neurons in net2')