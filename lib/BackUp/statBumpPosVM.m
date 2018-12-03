function [BumpPos, meanBumpPos, varBumpPos, concBumpPos] = statBumpPosVM(O, NetPars)
% Calculate the bump position by projecting network
% bump on an unit circle

% Wen-hao Zhang, Apr-30-2015
cirPos = exp(1i * NetPars.PrefStimDeg/ NetPars.WidthDeg * pi);

BumpPos = squeeze(mean(bsxfun(@times, cirPos, O), 1)); % average over all neurons
BumpPos = BumpPos ./ abs(BumpPos); % normalize
BumpPos(isnan(BumpPos)) = 0;

% Concentration of distribution of BumpPos
% For directional statistics, the concentration is the absolute value of summed points
meanBumpPos = mean(BumpPos, 2); % complex number

concBumpPos = abs(meanBumpPos); % mean resultant length, average over time
concBumpPos = mrl2Kappa(concBumpPos); % concentration parameter

% Angular
meanBumpPos = angle(meanBumpPos)* NetPars.WidthDeg/pi; % angular value
BumpPos = angle(BumpPos)* NetPars.WidthDeg/pi; % angular value

% Variance
devBumpPos = abs(bsxfun(@minus, BumpPos, meanBumpPos));
devBumpPos(devBumpPos>180) = devBumpPos(devBumpPos>180) - 360;
varBumpPos = sum(devBumpPos.^2, 2) / (size(devBumpPos,2)-1);


%%
% cirPos = exp(1i * NetPars.PrefStim/ NetPars.WidthDeg * pi);
%
% BumpPos = squeeze(sum(bsxfun(@times, cirPos, O), 1));
% BumpPos = angle(BumpPos)* NetPars.WidthDeg/pi;
