function [hFig, hAxe] = plotScanNetParsFunc(plotDatStruct, IdxPars, hFig)
% Plot the network performance of integration and keeping lost cue
% disparity inforamtion.

% Author: Wen-Hao Zhang, June-6-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University


%% Unfold the fields from plotDatStruct
for varName = fieldnames(plotDatStruct)';
    eval([varName{1} '=plotDatStruct.' varName{1} ';']);
end


%% Calculate the index of results to plot

if exist('IdxPars', 'var')
    for varName = fieldnames(IdxPars)';
        eval([varName{1} '=IdxPars.' varName{1} ';']);
    end
else
    % Without specified, Default is plotting all results
    IdxJrc = 1: length(NetPars.Jrc);
    IdxJrp = 1: length(NetPars.JrpRatio);
    IdxAmpl = 1: size(NetPars.Ampl, 2);
    IdxPosi = 1: size(NetPars.Posi, 2);
    IdxGroup = 1:4; % 1: congruent cell; 3: opposite cell
end


%% Generate Fig. and Axes handle
scrsz = get(0,'ScreenSize');
% close all;
if isempty('hFig')
    hFig = figure('Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
else
    set(hFig, 'Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
end

% figure(1);
for iter = 1: 6
    hAxe(iter) = subplot(2,3,iter);
    axis square;
    hold on
end

cSpec = colormap(jet(IdxAmpl(end)));
sSpec = 'o+o+';
% sSpec = 'o+sx';
lineSpec = {'-','-','--', '--'};
szMarker = 6;
lineWid = 1;

% Comparison between network results and Bayesian prediction
for iterGroup = IdxGroup
    for iterAmpl = IdxAmpl
        plot(hAxe(1), reshape(bayesResVM.meanNetBayes(iterAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(NetEstimRes.meanBumpPos(iterAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup),'color', cSpec(iterAmpl,:), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(2), reshape(bayesResVM.concNetBayes(iterAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            reshape(NetEstimRes.concBumpPos(iterAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup),'color', cSpec(iterAmpl,:), 'markersize', szMarker, 'linew', lineWid);
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
axis(180*[0 1 0 1])
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

%% Network results with position
IdxJrc = 1;
IdxJrp = 6;
IdxAmpl = 6;
IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell

hAxe(4) = subplot(4,3,7); hold on;
hAxe(5) = subplot(4,3,10); hold on;

disPos = diff(dimPar(4).valuePar);

for iterGroup = IdxGroup
    for iter = IdxAmpl
        % Concentration of congruent and opposite neurons
        plot(hAxe(4), disPos,...
            reshape(bayesResVM.concNetBayes(iter, IdxJrc, IdxJrp, :, 1,iterGroup),1,[]), ...
            'color', cSpec(iter,:), 'markerfacecolor', cSpec(iter,:), 'linew', lineWid);
        plot(hAxe(4), disPos,...
            reshape(NetEstimRes.concBumpPos(iter, IdxJrc, IdxJrp, :, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iter,:), 'linew', lineWid);
        
        % Mean of congruent and opposite neurons
        plot(hAxe(5), disPos,...
            reshape(bayesResVM.meanNetBayes(iter, IdxJrc, IdxJrp, :, 1,iterGroup),1,[]), ...
            'color', cSpec(iter,:), 'markerfacecolor', cSpec(iter,:), 'linew', lineWid);
        plot(hAxe(5), disPos,...
            reshape(NetEstimRes.meanBumpPos(iter, IdxJrc, IdxJrp, :, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iter,:), 'linew', lineWid);
        
        %  Firing rate with position
        plot(hAxe(3), disPos,...
            reshape(NetEstimRes.OHeightAvg(iter, IdxJrc, IdxJrp, :, 1, iterGroup),1,[]), ...
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
IdxJrc = 1;
IdxJrp = 1: 8;
IdxPosi = 2;

hAxe(6) = subplot(4,3,8); hold on;
hAxe(7) = subplot(4,3,11); hold on;

for iterGroup = IdxGroup
    for iter = IdxAmpl
        plot(hAxe(6), NetPars.JrpRatio(IdxJrp),...
            reshape(NetEstimRes.concBumpPos(IdxAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(6), NetPars.JrpRatio(IdxJrp),...
            reshape(bayesResVM.concNetBayes(IdxAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            'linew', lineWid);
        
        plot(hAxe(7), NetPars.JrpRatio(IdxJrp),...
            reshape(NetEstimRes.meanBumpPos(IdxAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
            sSpec(iterGroup), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(7), NetPars.JrpRatio(IdxJrp),...
            reshape(bayesResVM.meanNetBayes(IdxAmpl, IdxJrc, IdxJrp, IdxPosi, 1, iterGroup),1,[]), ...
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