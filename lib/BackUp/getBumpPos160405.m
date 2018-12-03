function [BumpPos, meanBumpPos, varBumpPos] = getBumpPos(O, NetPars)
% Calculate the bump position by projecting network
% bump on an unit circle

% Wen-hao Zhang, Apr-30-2015
cirPos = exp(1i * NetPars.PrefStim/ NetPars.Width * pi);

BumpPos = squeeze(sum(bsxfun(@times, cirPos, O), 1));
BumpPos = BumpPos ./ abs(BumpPos);
BumpPos(isnan(BumpPos)) = 0;

meanBumpPos = mean(BumpPos, 2);

% Angular
meanBumpPos = angle(meanBumpPos)* NetPars.Width/pi;
BumpPos = angle(BumpPos)* NetPars.Width/pi;

devBumpPos = abs(bsxfun(@minus, BumpPos, meanBumpPos));
devBumpPos(devBumpPos>180) = devBumpPos(devBumpPos>180) - 360;
varBumpPos = sum(devBumpPos.^2, 2) / (size(devBumpPos,2)-1);


% Concentration of parameter (von-Mises distribution)
% conBumpPos= 



%%
% cirPos = exp(1i * NetPars.PrefStim/ NetPars.Width * pi);
% 
% BumpPos = squeeze(sum(bsxfun(@times, cirPos, O), 1));
% BumpPos = angle(BumpPos)* NetPars.Width/pi;
