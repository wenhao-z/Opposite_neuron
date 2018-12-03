function NetStat = statNetBumpXTrial(O, NetPars)
% Calculate the statistics of network activities across TRIAL
% The statistics including
% 1. Bump position (by projecting network bump on an unit circle)
% 2. Bump height

% INPUT:
% O: [N, numGroup, Time, Trial]
%    The last dim of O must be TRIAL!!!

% Wen-hao Zhang, Nov-17-2016

%% Bump position
cirPos = exp(1i * NetPars.PrefStim/ NetPars.Width * pi);

% BumpPos = squeeze(mean(bsxfun(@times, cirPos, O), 1)); % average over all neurons
BumpPos = mean(bsxfun(@times, cirPos, O), 1); % average over all neurons
% BumpPos = shiftdim(BumpPos, 1);
BumpPos = BumpPos ./ abs(BumpPos); % normalize
BumpPos(isnan(BumpPos)) = 0;

% Concentration of distribution of BumpPos
% For directional statistics, the concentration is the absolute value of summed points
meanBumpPos = mean(BumpPos, 4); % complex number, average over the last dim (trial)

mrlBumpPos = abs(meanBumpPos); % mean resultant length, average over time
concBumpPos = mrl2Kappa(mrlBumpPos); % concentration parameter

% Angular
meanBumpPos = angle(meanBumpPos)* NetPars.Width/pi; % angular value
BumpPos = angle(BumpPos)* NetPars.Width/pi; % angular value

% Variance
devBumpPos = abs(bsxfun(@minus, BumpPos, meanBumpPos));
devBumpPos(devBumpPos>180) = devBumpPos(devBumpPos>180) - 360;
varBumpPos = sum(devBumpPos.^2, 4) / (size(devBumpPos,2)-1); % Average over the last dim(trial)

%% Bump height
OHeight = max(O, [], 1); % Average over neurons
OHeightAvgTrial = max( mean(O, 4), [], 1); % Average over neurons
% UHeightAvgTrial = max(U, [], 1); % Average over neurons


%% Fold statistics into a struct
NetStat.BumpPos = squeeze(BumpPos);
NetStat.meanBumpPos = squeeze(meanBumpPos);
NetStat.mrlBumpPos = squeeze(mrlBumpPos);
NetStat.concBumpPos = squeeze(concBumpPos);
NetStat.varBumpPos = squeeze(varBumpPos);

NetStat.OHeight = squeeze(OHeight);
NetStat.OHeightAvgTrial = squeeze(OHeightAvgTrial);

