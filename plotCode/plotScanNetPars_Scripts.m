% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;

setWorkPath;
datPath = fullfile(Path_RootDir, 'Data');
Folder = 'NetIntSeg';
% Folder = 'DoubleCANNs';

% fileName = 'scanNetPars_170606_2140.mat'; % krpRatio = 0.5;
% fileName = 'scanNetPars_170607_0310.mat'; % krpRatio = 0;
% fileName = 'scanNetPars_170607_1403.mat'; % krpRatio = 0.2;
% fileName = 'scanNetPars_170616_0200.mat';

% fileName = 'scanNetPars_170622_0049.mat';
% fileName = 'scanNetPars_170622_2134.mat';
% fileName = 'scanNetPars_170705_1819.mat'; % Noise reciprocal connections
% fileName = 'scanNetPars_170707_0646.mat';

% fileName = 'scanNetPars_180531_2234.mat'; % sameSeed, multi-trial
% fileName = 'scanNetPars_180604_2225.mat'; % sameSeed, single long trial
% fileName = 'scanNetPars_180605_1939.mat'; % sameSeed, multi-trial, new parameter range

% fileName = 'scanNetPars_180604_2215.mat'; % sameSeedCueCond, single long trial
fileName = 'scanNetPars_180613_2015.mat'; % sameSeedCueCond, multi-trial, new parameter range

plotDatStruct = load(fullfile(datPath, Folder, fileName), ...
    'NetEstimRes', 'bayesRes', 'dimPar', 'NetPars');


%% Permute bump position through adding random noise
bPermBumpPos = 1;
if bPermBumpPos
    szDat = size(plotDatStruct.bayesRes.meanNetBayes_VM);
    randArray = repmat(rand(szDat(1:end-1)), [ones(1, length(szDat)-1), szDat(end)]);
    randArray = (randArray - 0.5) * plotDatStruct.NetPars.Width*2;
    randArray = round(randArray);
    
    plotDatStruct.bayesRes.meanNetBayes_VM = plotDatStruct.bayesRes.meanNetBayes_VM + randArray;
    plotDatStruct.bayesRes.meanNetBayes_VM = ...
        angle(exp(1i * plotDatStruct.bayesRes.meanNetBayes_VM * pi/plotDatStruct.NetPars.Width))...
        * plotDatStruct.NetPars.Width/pi;
    
    namePar = {plotDatStruct.dimPar.namePar};
    IdxCueCond = find(cellfun(@(x) strcmp(x, 'cueCond'), namePar));
    szRep = ones(1, length(namePar)+1);
    szRep(IdxCueCond) = 3; % three cue conditions
    randArray = repmat(randArray, szRep);
    plotDatStruct.NetEstimRes.meanBumpPos = plotDatStruct.NetEstimRes.meanBumpPos + randArray;
    plotDatStruct.NetEstimRes.meanBumpPos = ...
        angle(exp(1i * plotDatStruct.NetEstimRes.meanBumpPos * pi/plotDatStruct.NetPars.Width))...
        * plotDatStruct.NetPars.Width/pi;
end
%% Parameters of plotted data
IdxPars_CmprDat.AmplRatio = 1:8;
IdxPars_CmprDat.JrcRatio = 1:2;
IdxPars_CmprDat.JrpRatio = 1:7; % 5;
IdxPars_CmprDat.IdxNeuronGroup = [3,4]; % 1: congruent cell; 3: opposite cell
IdxPars_CmprDat.Posi = 1:10;
IdxPars_CmprDat.krpRatio = 1:3;
IdxPars_CmprDat.cueCond = find(plotDatStruct.NetPars.cueCond == 0); % combined condition

dimSpecPlot.colorDimName = 'AmplRatio';
dimSpecPlot.colorDimName = 'IdxNeuronGroup';
dimSpecPlot.markerDimName = 'IdxNeuronGroup';

IdxPars_Posi.AmplRatio = 5; % multiple value allowed
IdxPars_Posi.JrcRatio = 1;
IdxPars_Posi.JrpRatio = 5;
IdxPars_Posi.IdxNeuronGroup = [1,3]; % 1: congruent cell; 3: opposite cell
IdxPars_Posi.krpRatio = 3;
IdxPars_Posi.cueCond = find(plotDatStruct.NetPars.cueCond == 0); % combined condition

IdxPars_JrpRatio.AmplRatio = 5; % multiple value allowed
IdxPars_JrpRatio.JrcRatio = 1;
IdxPars_JrpRatio.IdxNeuronGroup = [1,3]; % 1: congruent cell; 3: opposite cell
IdxPars_JrpRatio.krpRatio = 3;
IdxPars_JrpRatio.Posi = 9;
IdxPars_JrpRatio.cueCond = find(plotDatStruct.NetPars.cueCond == 0); % combined condition

IdxPars.CmprNetandBayes = IdxPars_CmprDat;
IdxPars.ResWithPosi = IdxPars_Posi;
IdxPars.ResWithJrp = IdxPars_JrpRatio;

[hFig, hAxe] = plotScanNetParsFunc(plotDatStruct, IdxPars, dimSpecPlot, figure(2));


%% Network results with input intensity
% fileName = 'scanNetPars_diffAmpl(12-Apr).mat';
%
% load(fullfile(datPath, Folder, fileName), ...
%     'meanNetSim', 'varNetSim', 'concNetSim', 'bayesRes', ...
%     'bayesResVM', 'OHeight', 'dimPar', 'NetPars');
%
% IdxJrc = 1;
% IdxJrp = 2;
%
% Ampl = unique(NetPars.Ampl)';
% IdxAmpl2 = 6;
% IdxAmpl = find(NetPars.Ampl(2,:)== Ampl(IdxAmpl2)); % change of cue 1 intensity
% IdxAmpl(2,:) = find(NetPars.Ampl(1,:)== Ampl(IdxAmpl2)); % change of cue 2 intensity
%
% IdxGroup = [1,3]; % 1: congruent cell; 3: opposite cell
%
% hAxe(8) = subplot(4,3,9); hold on;
% hAxe(9) = subplot(4,3,12); hold on;
%
% iterLayer = 1;
%
% for iterGroup = IdxGroup
%     for iter = IdxAmpl
%         % Concentration of congruent and opposite neurons
%         plot(hAxe(8), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
%             reshape(NetEstimRes.concBumpPos(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
%             sSpec(iterGroup), 'markersize', szMarker);
%         plot(hAxe(8), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
%             reshape(bayesResVM.concNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
%
%         plot(hAxe(9), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
%             reshape(NetEstimRes.meanBumpPos(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []), ...
%             sSpec(iterGroup), 'markersize', szMarker);
%         plot(hAxe(9), NetPars.Ampl(iterLayer, IdxAmpl(iterLayer,:))/NetPars.U0, ...
%             reshape(bayesResVM.meanNetBayes(IdxAmpl(iterLayer,:), IdxJrp,1, iterGroup), 1, []));
%
%     end
% end
% axes(hAxe(8));
% ylabel('Network concentration')
% axes(hAxe(9));
% xlabel('Intensity of cue 1/ (U_m^0)')
% ylabel('Network mean (\circ)')

%%
if bSavePlot
    cd([datPath, '/figure']);
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end