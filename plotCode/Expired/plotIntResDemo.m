% Plot the intermediate results
% Wen-Hao Zhang, 5-April-2016

hAxe = [];
cMax = max(InputSet.O(:));
tPlot = (1: size(InputSet.O, 3)) * NetPars.dt/NetPars.tau;
textTitle = {'Net 1, congruent', 'Net 2, congruent',...
    'Net 1, opposite', 'Net 2, opposite'};
for iter = 1: size(InputSet.O, 2)
    % Population activity across time
    hAxe = [hAxe, subplot(2, 8, (1:3) + (iter-1)*4)];
    imagesc(tPlot, NetPars.PrefStim, squeeze(InputSet.O(:,iter,:)))
    set(gca, 'ytick', [-178, -90: 90: 180], 'yticklabel', -180: 90: 180)
    ylabel('Neuron Index'); axis xy; 
%     colorbar
    caxis([0, cMax])
%     title
    title(textTitle{iter})
    
    % Temporal averaged population activity
    hAxe = [hAxe, subplot(2, 8, iter * 4)];
    plot(NetPars.PrefStim, OAvg(:,iter))
    hold on
    plot(netpars.Posi(1)*ones(1,2), [0, max(OAvg(:))], '--k')
    hold off
    axis tight; 
%     axis square
    set(gca, 'xtick', -180: 90: 180)
end
clear cMax

axes(hAxe(end-1))
colorbar


% Putative congruent neurons in network 1
% subplot(2,8,1:3)
% tPlot = (1: size(InputSet.O, 3)) * NetPars.dt/NetPars.tau;
% imagesc(tPlot, NetPars.PrefStim, squeeze(InputSet.O(:,1,:)))
% set(gca, 'ytick', [-178, -90: 90: 180], 'yticklabel', -180: 90: 180)
% ylabel('Neuron Index'); axis xy; colorbar
% 
% % Putative opposite neurons in network 1
% subplot(2,8,5:7)
% imagesc(tPlot, NetPars.PrefStim, squeeze(InputSet.O(:,2,:)))
% xlabel('Time (\tau)');
% ylabel('Neuron Index');axis xy; colorbar
% set(gca, 'ytick', [-178, -90: 90: 180], 'yticklabel', -180: 90: 180)
% 
% % Putative congruent neurons in network 1
% subplot(2,8, 4)
% plot(NetPars.PrefStim, OAvg(:,1))
% hold on
% plot(netpars.Posi(1)*ones(1,2), [min(OAvg(:,1)), max(OAvg(:,1))], '--k')
% hold off
% axis tight; axis square
% set(gca, 'xtick', -180: 90: 180)
% 
% % Putative opposite neurons in network 1
% subplot(2,8, 8)
% plot(NetPars.PrefStim, OAvg(:,2))
% hold on
% plot(netpars.Posi(2)*ones(1,2), [min(OAvg(:,2)), max(OAvg(:,2))], '--k')
% hold off
% axis tight; axis square
% xlabel('Neuron Index')
% set(gca, 'xtick', -180: 90: 180)