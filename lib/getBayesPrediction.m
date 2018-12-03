function bayesRes = getBayesPrediction(NetEstimRes, NetPars, dimPar)
% Calculate Bayesian predictions for integration and the lost disparity
% information by using von Mises distribution and Gaussian distribution

% Author: Wen-Hao Zhang, June-6-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

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

% Index to permute the array to move dim of cueCond into 1st dim.
% Permute the high-dim array into [cueCond, IdxGroup, ...]
% Calculate deviation needs to use IdxGroup information.
IdxPerm = ndims(NetEstimRes.meanBumpPos);
IdxPerm = [dimCueCond, IdxPerm, setdiff(1:IdxPerm, [dimCueCond, IdxPerm])];

%% Von Mises prediction
% Here, meanNetSim should be in rad!!!
meanNetSim = permute(NetEstimRes.meanBumpPos, IdxPerm)* pi/ NetPars.Width; % unit: rad
concNetSim = permute(NetEstimRes.concBumpPos, IdxPerm);

sizeArray = size(meanNetSim);
meanNetSim = reshape(meanNetSim, sizeArray(1), sizeArray(2), []);
concNetSim = reshape(concNetSim, sizeArray(1), sizeArray(2), []);

% Bayesian prediction
vecNetBayes = sum(concNetSim(end-1:end,:,:).* exp(1i* meanNetSim(end-1:end,:,:)), 1);

sizeArray(1) = 1; % two cue conditions --> 1 prediction, 
                           
% Reshape vecNetBayes as the same order and size of meanNetSim1
vecNetBayes = reshape(vecNetBayes, sizeArray);
vecNetBayes = ipermute(vecNetBayes, IdxPerm);

% Find the mean and concentration of vecNetBayes
meanNetBayes = angle(vecNetBayes) * NetPars.Width/pi;
concNetBayes = abs(vecNetBayes);

% Fold into an output struct
bayesRes.meanNetBayes_VM = meanNetBayes;
bayesRes.concNetBayes_VM = concNetBayes;

%% Gaussian prediction
meanNetSim = permute(NetEstimRes.meanBumpPos, IdxPerm); % unit: degree
varNetSim = permute(NetEstimRes.varBumpPos, IdxPerm);

sizeArray = size(meanNetSim);
meanNetSim = reshape(meanNetSim, sizeArray(1), sizeArray(2), []);
varNetSim = reshape(varNetSim, sizeArray(1), sizeArray(2), []);

% Bayesian prediction
varNetBayes = 1./sum(1./varNetSim(end-1:end,:,:),1);
weightInt = bsxfun(@rdivide, flip(varNetSim(end-1:end,:,:), 1), sum(varNetSim(end-1:end,:,:),1));
meanNetBayes = sum(weightInt .* meanNetSim(end-1:end,:,:), 1);

% Actual weight of direct cue of two networks (only works for two network modules)
wDirCueNetSim(:,[1,3],:)  = (meanNetSim(1,[1,3],:) - meanNetSim(3,[1,3],:))./(meanNetSim(2,[1,3],:) - meanNetSim(3,[1,3],:));
wDirCueNetSim(:,[2,4],:)  = (meanNetSim(1,[2,4],:) - meanNetSim(2,[2,4],:))./(meanNetSim(3,[2,4],:) - meanNetSim(2,[2,4],:));

% Predicted weights of direct cue of two networks
wDirCueNetBayes(:,[1,3],:) = varNetSim(3,[1,3],:)./ sum(varNetSim(end-1:end,[1,3],:),1);
wDirCueNetBayes(:,[2,4],:) = varNetSim(2,[2,4],:)./ sum(varNetSim(end-1:end,[2,4],:),1);

dMeanPct = wDirCueNetSim - wDirCueNetBayes;
dVarPerct = varNetSim(1,:,:) ./ varNetBayes - 1;

% Reshape
sizeArray(1) = 1;
bayesRes.meanNetBayes_Gauss    = ipermute(reshape(meanNetBayes, sizeArray), IdxPerm);
bayesRes.varNetBayes_Gauss     = ipermute(reshape(varNetBayes, sizeArray), IdxPerm);
bayesRes.wDirCueNetSim_Gauss   = ipermute(reshape(wDirCueNetSim, sizeArray), IdxPerm);
bayesRes.wDirCueNetBayes_Gauss = ipermute(reshape(wDirCueNetBayes, sizeArray), IdxPerm);
bayesRes.dMeanPct_Gauss        = ipermute(reshape(dMeanPct, sizeArray), IdxPerm);
bayesRes.dVarPerct_Gauss       = ipermute(reshape(dVarPerct, sizeArray), IdxPerm);

end