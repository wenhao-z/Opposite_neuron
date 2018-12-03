% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;


setWorkPath;
datPath = fullfile(Path_RootDir, 'Data');
Folder = 'NetIntSeg';
% fileName = 'scanNetPars_170623_0256.mat';
fileName = 'scanNetPars_180611_1334.mat';

load(fullfile(datPath, Folder, fileName));


%% Parameters of plotted data
IdxAmpl1 = 1:13;
IdxAmpl2 = 1:13;
[IdxAmpl1, IdxAmpl2] = meshgrid(IdxAmpl1, IdxAmpl2);
IdxAmpl = sub2ind(sqrt(size(NetPars.AmplRatio,2))*ones(1,2), IdxAmpl1, IdxAmpl2);
IdxAmpl = IdxAmpl(:)';
clear IdxAmpl1 IdxAmpl2

IdxPosi = 1:10;
IdxGroup = [1, 3]; % 1: congruent cell; 3: opposite cell


%% Generate Fig. and Axes handle
scrsz = get(0,'ScreenSize');
% close all;
% figure('Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6])

% figure(1);
for iter = 1: 6
    %     hAxe(iter) = subplot(3,4,iter);
    hAxe(iter) = subplot(2,3,iter);
    axis square;
    hold on
end

cSpec = colormap(jet(IdxAmpl(end)));
sSpec = 'oo++';
lineSpec = {'-','-','--', '--'};
szMarker = 6;
lineWid = 1;

%% Comparison between network results and Bayesian prediction
for iterGroup = IdxGroup
    for iter = IdxAmpl
        plot(hAxe(1), reshape(bayesRes.meanNetBayes_VM(iter, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(NetEstimRes.meanBumpPos(iter, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(2), reshape(bayesRes.concNetBayes_VM(iter, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(NetEstimRes.concBumpPos(iter, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
    end
end

axes(hAxe(1));
xlabel({'Bayesian mean', '(von-Mises)'})
ylabel('Network mean')
axis tight;
axisLim = 1.2* max(abs(axis));
plot(axisLim*[-1, 1], axisLim*[-1, 1], 'k--')
% plot(180*[0, 1], 180*[0, 1], 'k--')
axis tight;
% axis(180*[0 1 0 1])
set(hAxe(1), 'xtick', 0:45:180,...
    'ytick', 0:45:180)


axes(hAxe(2));
xlabel({'Bayesian concentration', '(von-Mises)'})
ylabel('Network concentration')

axis tight;
axisLim = [0.8* min(abs(axis)), 1.2* max(abs(axis))];
plot(axisLim, axisLim, 'k--')
axis tight;

set(hAxe(1:4), 'fontsize', 9)

%% Network results with input intensity
IdxPosi = 2;

Ampl = unique(NetPars.AmplRatio)';
IdxAmpl2 = 5;
IdxAmpl = find(NetPars.AmplRatio(2,:)== Ampl(IdxAmpl2)); % change of cue 1 intensity
IdxAmpl(2,:) = find(NetPars.AmplRatio(1,:)== Ampl(IdxAmpl2)); % change of cue 2 intensity

IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell

hAxe(3) = subplot(4,3,7); hold on;
hAxe(4) = subplot(4,3,10); hold on;

iterLayer = 1;

for iterGroup = IdxGroup
    %     for iter = IdxAmpl
    % Concentration of congruent and opposite neurons
    plot(hAxe(3), NetPars.AmplRatio(iterLayer, IdxAmpl(iterLayer,:)), ...
        reshape(NetEstimRes.concBumpPos(IdxAmpl(iterLayer,:), IdxPosi,1, iterGroup), 1, []), ...
        sSpec(iterGroup));
    plot(hAxe(3), NetPars.AmplRatio(iterLayer, IdxAmpl(iterLayer,:)), ...
        reshape(bayesRes.concNetBayes_VM(IdxAmpl(iterLayer,:), IdxPosi,1, iterGroup), 1, []));
    
    % Mean bump position
    plot(hAxe(4), NetPars.AmplRatio(iterLayer, IdxAmpl(iterLayer,:)), ...
        reshape(NetEstimRes.meanBumpPos(IdxAmpl(iterLayer,:), IdxPosi,1, iterGroup), 1, []), ...
        sSpec(iterGroup));
    plot(hAxe(4), NetPars.AmplRatio(iterLayer, IdxAmpl(iterLayer,:)), ...
        reshape(bayesRes.meanNetBayes_VM(IdxAmpl(iterLayer,:), IdxPosi,1, iterGroup), 1, []));
    
    % Bump height
    plot(hAxe(5), NetPars.AmplRatio(iterLayer, IdxAmpl(iterLayer,:)), ...
        reshape(NetEstimRes.OHeightAvg(IdxAmpl(iterLayer,:), IdxPosi,1, iterGroup), 1, []));
    
    %     end
end
axes(hAxe(3));
ylabel('Network concentration')
axis tight
xlim([0.25, 1.55])
set(gca, 'xtick', 0.3: 0.3:1.5)

textTitle = { ...
    ['AmplRatio2=' num2str(Ampl(IdxAmpl2)), ',' ...
    'JrcRatio=' num2str(NetPars.JrcRatio),] ...
    ['JrpRatio=' num2str(NetPars.JrpRatio), ',' ...
    'krpRatio=' num2str(NetPars.krpRatio)]};
title(textTitle)

axes(hAxe(4));
xlabel('Intensity of cue / (U_m^0)')
ylabel('Network mean')
axis tight
ylim([-20, 20])
xlim([0.25, 1.55])
set(gca, 'xtick', 0.3: 0.3:1.5)

axes(hAxe(5))
xlabel('Intensity of cue / (U_m^0)')
ylabel('Bump height (Hz)')
xlim([0.25, 1.55])
ylim([0, 60])
set(gca, 'xtick', 0.3: 0.3:1.5)
title(textTitle)

%%
if bSavePlot
    cd([datPath, '/figure']);
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end