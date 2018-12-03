% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;

datPath = 'D:\Projects\InfoSegregation\DecenNet';
fileName = 'scanNetPars(07-Apr).mat';

load(fullfile(datPath, fileName));

%% Parameters of plots

% Name of all scanned parameters
namePar = {dimPar.namePar};

% Index of conditioned parameters but which cannot be plotted on the figure
IdxStruct.JrpRatio = 2;
IdxStruct.Ampl = 2;
IdxStruct.Posi = 3;
IdxStruct.cueCond = 1;

% mapping above index as the order of meanNetSim
IdxSet = size(meanNetSim); % Index for plotted data
IdxSet = mat2cell(IdxSet, 1, ones(1, length(IdxSet)) );
nameStruct = fieldnames(IdxStruct);
for iter = 1: length(nameStruct)
    Idx = cellfun(@(x) strcmp(x, namePar(iter)), nameStruct);
    Idx = find(Idx == 1);
    IdxSet{iter} = IdxStruct.(nameStruct{Idx});
end

IdxStruct.Group = [1,3];
IdxSet{end} = IdxStruct.Group;

% Index of parameters shown on the figure
% Parameter of x-axis
xPar = 'Posi';
xParDim = cellfun(@(x) strcmp(x, xPar), namePar);
xParDim = find(xParDim == 1);

% Parameter indicated by COLOR of dots or lines
clPar = 'Ampl';
clParDim = cellfun(@(x) strcmp(x, clPar), namePar);
clParDim = find(clParDim == 1);

IdxSet{xParDim} = 1:size(dimPar(xParDim).valuePar, 2);
IdxSet{clParDim} = 1:size(dimPar(clParDim).valuePar, 2);

%% Use IdxSet to construct plotStruct

expr = [];
for iter = 1: length(IdxSet)
    expr = [expr, '[IdxSet{' num2str(iter) '}],'];
end
expr = ['(' expr(1:end-1) ')'];
plotStruct.meanNetSim = squeeze(eval(['meanNetSim' expr])  );
plotStruct.varNetSim = squeeze(eval(['varNetSim' expr])  );
plotStruct.concNetSim = squeeze(eval(['concNetSim' expr])  );

if xParDim > clParDim
    plotStruct.meanNetSim = permute(plotStruct.meanNetSim, [2,1,3]);
    plotStruct.varNetSim = permute(plotStruct.varNetSim, [2,1,3]);
    plotStruct.concNetSim = permute(plotStruct.concNetSim, [2,1,3]);
end

%% Plot

for iter = xPar

end

%%
if bSavePlot
    cd([datPath, '/figure']); 
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end