function bayesRes = getBayesPredictionVM(meanNetSim, concNetSim, NetPars, dimPar)
% Calculate Bayesian predictions for integrated estimation when treating
% the samples are from von-Mises distribution
% Wen-Hao Zhang, 5-April-2016

% NOTE: meanNetSim and concNetSim must be in RADIAN!!!
% This code could be only used for two coupled nets!!
% The order of group: 1c, 2c, 1o, 2o

%% Bayesian estimation of mean and variance

for iter = 1: length(dimPar)
    if strcmp(dimPar(iter).namePar, 'cueCond')
        dimCueCond = iter;
        break
    end
end
clear iter

%% Permute the high-dim array into [Net group Idx, CudCond, ...]
IdxPerm = ndims(meanNetSim);
IdxPerm = [IdxPerm, dimCueCond, setdiff(1:IdxPerm, [IdxPerm, dimCueCond])];

meanNetSim1 = permute(meanNetSim, IdxPerm);
concNetSim1 = permute(concNetSim, IdxPerm);

sizeArray = size(meanNetSim1);
meanNetSim1 = reshape(meanNetSim1, sizeArray(1), sizeArray(2), []);
concNetSim1 = reshape(concNetSim1, sizeArray(1), sizeArray(2), []);

% Throw away the statistics under combined cue condition
meanNetSim1 = meanNetSim1(:, end-1:end, :);
concNetSim1 = concNetSim1(:, end-1:end, :);

sizeArray(2) = 1; % two cue conditions --> 1 prediction, 
                           % left for later use
%% Bayesian prediction
vecNetBayes = sum(concNetSim1 .* exp(1i* meanNetSim1), 2);

% Reshape vecNetBayes as the same order and size of meanNetSim1
vecNetBayes = reshape(vecNetBayes, sizeArray);
vecNetBayes = ipermute(vecNetBayes, IdxPerm);

% Find the mean and concentration of vecNetBayes
meanNetBayes = angle(vecNetBayes) * NetPars.Width/pi;
concNetBayes = abs(vecNetBayes);

%% Fold into an output struct
bayesRes.concNetBayes = concNetBayes;
bayesRes.meanNetBayes = meanNetBayes;


end