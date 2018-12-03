% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;

setWorkPath;
datPath = fullfile(Path_RootDir, 'Data');
Folder = 'Gauss';
fileName = 'scanNetPars_diffAmpl(12-Apr).mat';

load(fullfile(datPath, Folder, fileName), ...
    'meanNetSim', 'varNetSim', 'concNetSim', 'bayesRes', ...
    'bayesResVM', 'OHeight', 'dimPar', 'NetPars');


%% Parameters of plotted data
IdxJrc = 1;
IdxJrp = 1:3;

% IdxAmpl = 1: size(NetPars.Ampl,2);

IdxAmpl1 = 8:21;
IdxAmpl2 = 8:21;
[IdxAmpl1, IdxAmpl2] = meshgrid(IdxAmpl1, IdxAmpl2);
IdxAmpl = sub2ind(sqrt(size(NetPars.Ampl,2))*ones(1,2), IdxAmpl1, IdxAmpl2);
IdxAmpl = IdxAmpl(:)';
clear IdxAmpl1 IdxAmpl2

% IdxGroup = 1:4; % 1: congruent cell; 3: opposite cell
IdxGroup = [1, 3]; % 1: congruent cell; 3: opposite cell
% IdxGroup = [2, 4]; % 1: congruent cell; 3: opposite cell


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
        plot(hAxe(1), reshape(bayesResVM.meanNetBayes(iter, IdxJrp, 1, iterGroup),1,[]), ...
            reshape(meanNetSim(iter, IdxJrp, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(2), reshape(bayesResVM.concNetBayes(iter, IdxJrp, 1, iterGroup),1,[]), ...
            reshape(concNetSim(iter, IdxJrp, 1, iterGroup),1,[]), ...
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
IdxJrp = 2;

Ampl = unique(NetPars.Ampl)';
IdxAmpl2 = 6;
IdxAmpl = find(NetPars.Ampl(2,:)== Ampl(IdxAmpl2)); % change of cue 1 intensity
IdxAmpl(2,:) = find(NetPars.Ampl(1,:)== Ampl(IdxAmpl2)); % change of cue 2 intensity

IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell

hAxe(3) = subplot(4,3,7); hold on;
hAxe(4) = subplot(4,3,10); hold on;

iterLayer = 1;

for iterGroup = IdxGroup
    for iter = IdxAmpl
        % Concentration of congruent and opposite neurons        
        plot(hAxe(3), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(concNetSim(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
            sSpec(iterGroup));
        plot(hAxe(3), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(bayesResVM.concNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
        
        
        plot(hAxe(4), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(meanNetSim(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
            sSpec(iterGroup));
        plot(hAxe(4), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
            reshape(bayesResVM.meanNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
        
    end
end
axes(hAxe(3));
ylabel('Network concentration')
axes(hAxe(4));
xlabel('Intensity of cue / (U_m^0)')
ylabel('Network mean')

%%
if bSavePlot
    cd([datPath, '/figure']); 
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end