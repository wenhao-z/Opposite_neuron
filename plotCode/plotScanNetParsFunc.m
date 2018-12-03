function [hFig, hAxe] = plotScanNetParsFunc(plotDatStruct, IdxPars, dimSpecPlot, hFig)
% Plot the network performance of integration and keeping lost cue
% disparity inforamtion.

% Author: Wen-Hao Zhang, June-6-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University


%% Unfold the fields from struct
for varName = fieldnames(plotDatStruct)';
    eval([varName{1} '=plotDatStruct.' varName{1} ';']);
end

IdxPars_CmprDat = IdxPars.CmprNetandBayes;
IdxPars_Posi = IdxPars.ResWithPosi;
IdxPars_JrpRatio = IdxPars.ResWithJrp;

%% Generate Fig. and Axes handle
scrsz = get(0,'ScreenSize');
% close all;
if isempty('hFig')
    hFig = figure('Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
else
    set(hFig, 'Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
    clf;
end

for iter = 1: 6
    hAxe(iter) = subplot(2,3,iter);
    axis square;
    hold on
end

cSpec = colormap(jet(length(IdxPars_CmprDat.(dimSpecPlot.colorDimName))));
sSpec = 'o+o+';
% sSpec = 'o+sx';
lineSpec = {'-','-','--', '--'};
szMarker = 6;
lineWid = 1;

%% Crop the high-dim array
% Get the meaning of each dim of parGrid
namePar = {dimPar.namePar};
nameTestRes_Crop = {'bayesRes.meanNetBayes_VM', ...
    'bayesRes.concNetBayes_VM', ...
    'NetEstimRes.meanBumpPos', ...
    'NetEstimRes.concBumpPos', ...
    'NetEstimRes.OHeightAvg'};

%% Comparison between network results and Bayesian prediction
if ~exist('IdxPars_CmprDat', 'var')
    for varName = fieldnames(dimPar)
        IdxPars_CmprDat.(varName{1}) = 1: size(NetPars.(varName{1}), 2);
    end
    IdxPars_CmprDat.IdxNeuronGroup = 1:4; % 1: congruent cell; 3: opposite cell
    IdxPars_CmprDat.cueCond = find(NetPars.cueCond == 0); % Combined condition
end

% The index of dim in parGrid demonstrated in color and marker
if strcmp(dimSpecPlot.colorDimName, 'IdxNeuronGroup')
    % Caution: make sure the LAST dim is index of group (check with the code scanNetPars.)
    IdxColorDim = length(namePar) + 1;
else
    IdxColorDim = find(cellfun(@(x) strcmp(x, dimSpecPlot.colorDimName), namePar));
end
if strcmp(dimSpecPlot.markerDimName, 'IdxNeuronGroup')
    % Caution: make sure the LAST dim is index of group (check with the code scanNetPars.)
    IdxMarkerDim = length(namePar) + 1;
else
    IdxMarkerDim = find(cellfun(@(x) strcmp(x, dimSpecPlot.markerDimName), namePar));
end

if IdxMarkerDim == IdxColorDim
    for iter = 1: length(nameTestRes_Crop)
        TestRes = eval(nameTestRes_Crop{iter});
        TestRes = cropHighDimArray(TestRes, namePar, IdxPars_CmprDat, IdxMarkerDim);
        eval([nameTestRes_Crop{iter} '= TestRes;']);
    end
    Idx = ((bayesRes.meanNetBayes_VM(:) - NetEstimRes.meanBumpPos(:)) > NetPars.Width);
    bayesRes.meanNetBayes_VM(Idx) = bayesRes.meanNetBayes_VM(Idx) - 2*NetPars.Width;
    Idx = ((bayesRes.meanNetBayes_VM(:) - NetEstimRes.meanBumpPos(:)) < -NetPars.Width);
    bayesRes.meanNetBayes_VM(Idx) = bayesRes.meanNetBayes_VM(Idx) + 2*NetPars.Width;
    
    for iterMarkerDim = 1: length(IdxPars_CmprDat.(dimSpecPlot.markerDimName))
        plot(hAxe(1), squeeze(bayesRes.meanNetBayes_VM(iterMarkerDim, :)), ...
            squeeze(NetEstimRes.meanBumpPos(iterMarkerDim, :)), ...
            sSpec(iterMarkerDim),'color', cSpec(iterMarkerDim,:), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(2), squeeze(bayesRes.concNetBayes_VM(iterMarkerDim, :)), ...
            squeeze(NetEstimRes.concBumpPos(iterMarkerDim, :)), ...
            sSpec(iterMarkerDim),'color', cSpec(iterMarkerDim,:), 'markersize', szMarker, 'linew', lineWid);
    end
    
else
    for iter = 1: length(nameTestRes_Crop)
        TestRes = eval(nameTestRes_Crop{iter});
        TestRes = cropHighDimArray(TestRes, namePar, IdxPars_CmprDat, [IdxMarkerDim, IdxColorDim]);
        eval([nameTestRes_Crop{iter} '= TestRes;']);
    end
    
    for iterMarkerDim = 1: length(IdxPars_CmprDat.(dimSpecPlot.markerDimName))
        for iterColorDim = 1: length(IdxPars_CmprDat.(dimSpecPlot.colorDimName))
            plot(hAxe(1), squeeze(bayesRes.meanNetBayes_VM(iterMarkerDim, iterColorDim, :)), ...
                squeeze(NetEstimRes.meanBumpPos(iterMarkerDim, iterColorDim, :)), ...
                sSpec(iterMarkerDim),'color', cSpec(iterColorDim,:), 'markersize', szMarker, 'linew', lineWid);
            plot(hAxe(2), squeeze(bayesRes.concNetBayes_VM(iterMarkerDim, iterColorDim, :)), ...
                squeeze(NetEstimRes.concBumpPos(iterMarkerDim, iterColorDim, :)), ...
                sSpec(iterMarkerDim),'color', cSpec(iterColorDim,:), 'markersize', szMarker, 'linew', lineWid);
        end
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
set(hAxe(1), 'xtick', 0:45:180,'ytick', 0:45:180)
% Estimate the R-square to evaluate the accordance between network
% and Bayesian prediction.
[~, gof] = fit(bayesRes.meanNetBayes_VM(:), NetEstimRes.meanBumpPos(:), 'poly1');
axisLim = axis;
text(axisLim(2), axisLim(3), ['R^2=', num2str(gof.rsquare)], 'horizontalalignment', 'right', ...
    'verticalalignment', 'bottom');
title(...
{['AmplRatio=[', num2str(NetPars.AmplRatio(1,IdxPars_CmprDat.AmplRatio([1,end]))) , ...
'], JrcRatio=[' num2str(NetPars.JrcRatio(IdxPars_CmprDat.JrcRatio([1,end]))), ']' ], ...
['JrpRatio=[' num2str(NetPars.JrpRatio(IdxPars_CmprDat.JrpRatio([1,end]))), ...
'], krpRatio=[' num2str(NetPars.krpRatio(IdxPars_CmprDat.krpRatio([1,end]))), ']'  ]});


axes(hAxe(2));
xlabel({'Bayesian concentration', '(von-Mises)'})
ylabel('Network concentration')
axis tight;
axisLim = [0.8* min(abs(axis)), 1.2* max(abs(axis))];
plot(axisLim, axisLim, 'k--')
axis tight;
% Estimate the R-square to evaluate the accordance between network
% and Bayesian prediction.
[~, gof] = fit(bayesRes.concNetBayes_VM(:), NetEstimRes.concBumpPos(:), 'poly1');
axisLim = axis;
text(axisLim(2), axisLim(3), ['R^2=', num2str(gof.rsquare)], 'horizontalalignment', 'right', ...
    'verticalalignment', 'bottom');
clear gof

set(hAxe(1:4), 'fontsize', 9)
%% Network results with position
% ================================================================================================
hAxe(4) = subplot(4,3,7); hold on;
hAxe(5) = subplot(4,3,10); hold on;

disPos = diff(NetPars.Posi(1:2,:), 1, 1);

IdxPosiDim = find(cellfun(@(x) strcmp(x, 'Posi'), namePar));
IdxAmplDim = find(cellfun(@(x) strcmp(x, 'AmplRatio'), namePar));
IdxNeuronGroup = length(namePar) + 1;

% Unfold the fields from plotDatStruct again
for varName = fieldnames(plotDatStruct)';
    eval([varName{1} '=plotDatStruct.' varName{1} ';']);
end
for iter = 1: length(nameTestRes_Crop)
    TestRes = eval(nameTestRes_Crop{iter});
    % Crop the high-dim array on each dim of parameters
    TestRes = cropHighDimArray(TestRes, namePar, IdxPars_Posi, [IdxPosiDim, IdxAmplDim, IdxNeuronGroup]);
    eval([nameTestRes_Crop{iter} '= TestRes;']);
end

for iterGroup = 1: length(IdxPars_Posi.IdxNeuronGroup)
    for iterAmpl = 1: length(IdxPars_Posi.AmplRatio)
        % Concentration of congruent and opposite neurons
        plot(hAxe(4), disPos,...
            bayesRes.concNetBayes_VM(:,iterAmpl, iterGroup), ...
            'color', cSpec(iterAmpl,:), 'markerfacecolor', cSpec(iterAmpl,:), 'linew', lineWid);
        plot(hAxe(4), disPos,...
            NetEstimRes.concBumpPos(:,iterAmpl, iterGroup), ...
            sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iterAmpl,:), 'linew', lineWid);
        
        % Mean of congruent and opposite neurons
        plot(hAxe(5), disPos,...
            bayesRes.meanNetBayes_VM(:,iterAmpl, iterGroup), ...
            'color', cSpec(iterAmpl,:), 'markerfacecolor', cSpec(iterAmpl,:), 'linew', lineWid);
        plot(hAxe(5), disPos,...
            NetEstimRes.meanBumpPos(:,iterAmpl, iterGroup), ...
            sSpec(iterGroup), 'markersize', szMarker, 'color', cSpec(iterAmpl,:), 'linew', lineWid);
        
        %  Firing rate with position
        plot(hAxe(3), disPos,...
            NetEstimRes.OHeightAvg(:,iterAmpl, iterGroup), ...
            'linestyle', lineSpec{iterGroup}, 'color', cSpec(iterAmpl,:), 'linew', lineWid);
    end
end

textTitle = { ...
['AmplRatio=', num2str(NetPars.AmplRatio(1,IdxPars_Posi.AmplRatio)), ', ' ...
'JrcRatio=' num2str(NetPars.JrcRatio(IdxPars_Posi.JrcRatio)), ', '] ...
['JrpRatio=' num2str(NetPars.JrpRatio(IdxPars_Posi.JrpRatio)), ', ' ...
'krpRatio=' num2str(NetPars.krpRatio(IdxPars_Posi.krpRatio))]};

axes(hAxe(4));
% xlabel({'Cue disparity', 'x_2 - x_1'})
ylabel('Network concentration')
title(textTitle)
axes(hAxe(5));
xlabel({'Cue disparity', 'x_2 - x_1'})
ylabel('Network mean (\circ)')
axes(hAxe(3))
ylabel('Firing rate (Hz)')
xlabel({'Cue disparity', 'x_2 - x_1'})
title(textTitle)
set(hAxe(3:5), 'xlim', [0, 180], 'xtick', 0:45:180)

%%  Net results with reciprocal connection
hAxe(6) = subplot(4,3,8); hold on;
hAxe(7) = subplot(4,3,11); hold on;

IdxJrpDim = find(cellfun(@(x) strcmp(x, 'JrpRatio'), namePar));

% Unfold the fields from plotDatStruct again
for varName = fieldnames(plotDatStruct)';
    eval([varName{1} '=plotDatStruct.' varName{1} ';']);
end
for iter = 1: length(nameTestRes_Crop)
    TestRes = eval(nameTestRes_Crop{iter});
    % Crop the high-dim array on each dim of parameters
    TestRes = cropHighDimArray(TestRes, namePar, IdxPars_JrpRatio, [IdxJrpDim, IdxAmplDim, IdxNeuronGroup]);
    eval([nameTestRes_Crop{iter} '= TestRes;']);
end

for iterGroup = 1: length(IdxPars_Posi.IdxNeuronGroup)
    for iterAmpl = 1: length(IdxPars_Posi.AmplRatio)
        plot(hAxe(6), NetPars.JrpRatio,...
            NetEstimRes.concBumpPos(:,iterAmpl, iterGroup), ...
            sSpec(iterGroup), 'color', cSpec(iterAmpl,:),'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(6), NetPars.JrpRatio,...
            bayesRes.concNetBayes_VM(:,iterAmpl, iterGroup), ...
            'color', cSpec(iterAmpl,:), 'linew', lineWid);
        
        plot(hAxe(7), NetPars.JrpRatio,...
            NetEstimRes.meanBumpPos(:,iterAmpl, iterGroup), ...
            sSpec(iterGroup), 'color', cSpec(iterAmpl,:), 'markersize', szMarker, 'linew', lineWid);
        plot(hAxe(7), NetPars.JrpRatio,...
            bayesRes.meanNetBayes_VM(:,iterAmpl, iterGroup), ...
            'color', cSpec(iterAmpl,:), 'linew', lineWid);
    end
end

textTitle = { ...
['AmplRatio=', num2str(NetPars.AmplRatio(1,IdxPars_JrpRatio.AmplRatio)), ', ' ...
'JrcRatio=' num2str(NetPars.JrcRatio(IdxPars_JrpRatio.JrcRatio)), ', '] ...
['diffPosi=' num2str(diff(NetPars.Posi(1:2,IdxPars_JrpRatio.Posi))), ', ' ...
'krpRatio=' num2str(NetPars.krpRatio(IdxPars_JrpRatio.krpRatio))]};

axes(hAxe(6))
ylabel('Network concentration')
title(textTitle)
axes(hAxe(7))
xlabel('Reciprocal strength J_{rp} / (J_{rc})')
ylabel('Network mean (\circ)')
set(hAxe(6:7), 'xlim', NetPars.JrpRatio([1, end]), ...
    'xtick', 0.3: 0.3: 0.9)

end