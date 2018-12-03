function [NetStat, BumpPos] = statNetBumpXTime(O, NetPars)
% Calculate the bump position by projecting network
% bump on an unit circle

% INPUT
% O: [N, numGroup, Time, Trial]
% The statistics of bump positon is averaged over time.

% Wen-hao Zhang, Nov-17-2016

cirPos = exp(1i * NetPars.PrefStim/ NetPars.Width * pi);

BumpPos = mean(bsxfun(@times, cirPos, O), 1); % average over all neurons
BumpPos = shiftdim(BumpPos, 1);
BumpPos = BumpPos ./ abs(BumpPos); % normalize
BumpPos(isnan(BumpPos)) = 0;

% Concentration of distribution of BumpPos
% For directional statistics, the concentration is the absolute value of summed points
meanBumpPos = mean(BumpPos, 2); % complex number

mrlBumpPos = abs(meanBumpPos); % mean resultant length, average over time
concBumpPos = mrl2Kappa(mrlBumpPos); % concentration parameter

% Angular
meanBumpPos = angle(meanBumpPos)* NetPars.Width/pi; % angular value
BumpPos = angle(BumpPos)* NetPars.Width/pi; % angular value

% Variance
devBumpPos = abs(bsxfun(@minus, BumpPos, meanBumpPos));
devBumpPos(devBumpPos>180) = devBumpPos(devBumpPos>180) - 360;
varBumpPos = sum(devBumpPos.^2, 2) / (size(devBumpPos,2)-1);

% Bump height
OHeight = max(mean(O, 3), [], 1); % Average over neurons
% UHeight = max(mean(U, 3), [], 1);

%% Fold the results into a struct
% NetStat.BumpPos = squeeze(BumpPos);
NetStat.meanBumpPos = squeeze(meanBumpPos);
NetStat.mrlBumpPos = squeeze(mrlBumpPos);
NetStat.concBumpPos = squeeze(concBumpPos);
NetStat.varBumpPos = squeeze(varBumpPos);

NetStat.OHeight = squeeze(OHeight)';

