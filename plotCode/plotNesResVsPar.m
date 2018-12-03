% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;

datPath = 'D:\Projects\InfoSegregation\DecenNet\Data\';
Folder = 'VMTrial';
% Folder = [];
% fileName = 'scanNetPars_Jrc4e-1.mat';
% fileName = 'scanNetPars(11-Apr).mat';
% fileName = 'scanNetPars(10-May)_Jrc4e-1.mat';
fileName = 'scanNetPars(10-May).mat';

load(fullfile(datPath, Folder, fileName), ...
    'meanNetSim', 'varNetSim', 'concNetSim', 'bayesRes', ...
    'bayesResVM', 'OHeight', 'dimPar', 'NetPars');

%% Parameters of plotted data
IdxJrc = 1;
IdxJrp = 1:9;
IdxAmpl = 1:12;
IdxPosi = 1:10;
IdxGroup = 1:4; % 1: congruent cell; 3: opposite cell
IdxGroup = [1, 3]; % 1: congruent cell; 3: opposite cell
% IdxGroup = [2, 4]; % 1: congruent cell; 3: opposite cell
% IdxGroup = [1,2]; % 1: congruent cell; 3: opposite cell

%% Generate Fig. and Axes handle
scrsz = get(0,'ScreenSize');
% close all;
figure('Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6])

% figure(1);
for iter = 1: 6
%     hAxe(iter) = subplot(3,4,iter);
    hAxe(iter) = subplot(2,3,iter);
    axis square;
    hold on
end

cSpec = colormap(jet(IdxAmpl(end)));
% sSpec = 'oo++';
sSpec = 'os+x';
lineSpec = {'-','-','--', '--'};
szMarker = 6;
lineWid = 1;

%% Comparison between network results and Bayesian prediction
for iterGroup = IdxGroup
    for iter = IdxAmpl
        plot(hAxe(1), reshape(bayesResVM.meanNetBayes(iter, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(meanNetSim(iter, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup),'color', cSpec(iter,:), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(2), reshape(bayesResVM.concNetBayes(iter, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(concNetSim(iter, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup),'color', cSpec(iter,:), 'markersize', szMarker, 'linew', lineWid);
    end
end

axes(hAxe(1)); 
xlabel({'Bayesian mean (\circ)', '(von-Mises)'})
ylabel('Network mean (\circ)')
axis tight;
axisLim = 1.2* max(abs(axis));
plot(axisLim*[-1, 1], axisLim*[-1, 1], 'k--')
% plot(180*[0, 1], 180*[0, 1], 'k--')
axis tight;
% axis(180*[0 1 0 1])
% set(hAxe(1), 'xtick', 0:45:180,...
%     'ytick', 0:45:180)


axes(hAxe(2)); 
xlabel({'Bayesian concentration', '(von-Mises)'})
ylabel('Network concentration')
axis tight;
axisLim = [0.8* min(abs(axis)), 1.2* max(abs(axis))];
plot(axisLim, axisLim, 'k--')
axis tight;
% set(hAxe(2), 'xscale', 'log', 'yscale', 'log')
set(hAxe(1:4), 'fontsize', 9)

%% Network results with position
IdxJrp = 7;
IdxAmpl = 12;
IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell

hAxe(4) = subplot(4,3,7); hold on;
hAxe(5) = subplot(4,3,10); hold on;

disPos = diff(dimPar(3).valuePar);

for iterGroup = IdxGroup
    for iter = IdxAmpl
        % Concentration of congruent and opposite neurons
        plot(hAxe(4), disPos,...
            reshape(bayesResVM.concNetBayes(iter, IdxJrp, :, 1,iterGroup),1,[]), ...
            'color', cSpec(iter,:), 'markerfacecolor', cSpec(iter,:), 'linew', lineWid);
        plot(hAxe(4), disPos,...
            reshape(concNetSim(iter, IdxJrp, :, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iter,:), 'linew', lineWid);
        
        % Mean of congruent and opposite neurons
        plot(hAxe(5), disPos,...
            reshape(bayesResVM.meanNetBayes(iter, IdxJrp, :, 1,iterGroup),1,[]), ...
            'color', cSpec(iter,:), 'markerfacecolor', cSpec(iter,:), 'linew', lineWid);
        plot(hAxe(5), disPos,...
            reshape(meanNetSim(iter, IdxJrp, :, 1, iterGroup),1,[]), ...
           sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iter,:), 'linew', lineWid);
        
        %  Firing rate with position
        plot(hAxe(3), disPos,...
            reshape(OHeight(iter, IdxJrp, :, 1, iterGroup),1,[]), ...
            'linestyle', lineSpec{iterGroup}, 'color', cSpec(iter,:), 'linew', lineWid);
    end
end
axes(hAxe(4));
% xlabel({'Cue disparity', 'x_2 - x_1'})
ylabel('Network concentration')
axes(hAxe(5));
xlabel({'Cue disparity', 'x_2 - x_1'})
ylabel('Network mean (\circ)')
axes(hAxe(3))
ylabel('Firing rate (Hz)')
xlabel({'Cue disparity', 'x_2 - x_1'})

set(hAxe(3:5), 'xlim', [0, 180], 'xtick', 0:45:180)

%%  Net results with reciprocal connection
IdxJrp = 1: 8;
IdxPosi = 3;

hAxe(6) = subplot(4,3,8); hold on;
hAxe(7) = subplot(4,3,11); hold on;

for iterGroup = IdxGroup
    for iter = IdxAmpl
        plot(hAxe(6), NetPars.JrpRatio(IdxJrp),...
            reshape(concNetSim(IdxAmpl, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(6), NetPars.JrpRatio(IdxJrp),...
            reshape(bayesResVM.concNetBayes(IdxAmpl, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            'linew', lineWid);
       
        plot(hAxe(7), NetPars.JrpRatio(IdxJrp),...
            reshape(meanNetSim(IdxAmpl, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(7), NetPars.JrpRatio(IdxJrp),...
            reshape(bayesResVM.meanNetBayes(IdxAmpl, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            'linew', lineWid);
    end
end
axes(hAxe(6))
% xlabel('Reciprocal strength J_{rp} / J_{rc}')
ylabel('Network concentration')
axes(hAxe(7))
xlabel('Reciprocal strength J_{rp} / (J_{rc})')
ylabel('Network mean (\circ)')
set(hAxe(6:7), 'xlim', NetPars.JrpRatio([IdxJrp(1), IdxJrp(end)]), ...
    'xtick', 0.3: 0.3: 0.9)

%% Network results with input intensity
fileName = 'scanNetPars_diffAmpl(12-Apr).mat';
Folder = 'Gauss';
load(fullfile(datPath, Folder, fileName), ...
    'meanNetSim', 'varNetSim', 'concNetSim', 'bayesRes', ...
    'bayesResVM', 'OHeight', 'dimPar', 'NetPars');

IdxJrp = 2;

Ampl = unique(NetPars.Ampl)';
IdxAmpl2 = 6;
IdxAmpl = find(NetPars.Ampl(2,:)== Ampl(IdxAmpl2)); % change of cue 1 intensity
IdxAmpl(2,:) = find(NetPars.Ampl(1,:)== Ampl(IdxAmpl2)); % change of cue 2 intensity

IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell

hAxe(8) = subplot(4,3,9); hold on;
hAxe(9) = subplot(4,3,12); hold on;

iterLayer = 1;

for iterGroup = IdxGroup
    for iter = IdxAmpl
        % Concentration of congruent and opposite neurons        
        plot(hAxe(8), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(concNetSim(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
            sSpec(iterGroup), 'markersize', szMarker);
        plot(hAxe(8), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(bayesResVM.concNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
        
        plot(hAxe(9), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(meanNetSim(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
            sSpec(iterGroup), 'markersize', szMarker);
        plot(hAxe(9), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(bayesResVM.meanNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
        
    end
end
axes(hAxe(8));
ylabel('Network concentration')
axes(hAxe(9));
xlabel('Intensity of cue 1/ (U_m^0)')
ylabel('Network mean (\circ)')

%%
if bSavePlot
    cd([datPath, '/figure']); 
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end