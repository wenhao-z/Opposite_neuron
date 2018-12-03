function [BumpPos, meanBumpPos, varBumpPos, concBumpPos, mrlBumpPos] = statBumpPos(Input, NetPars, inputName)
% Calculate the bump position by projecting network
% bump on an unit circle

% Wen-hao Zhang, Dec-31-2016

if nargin < 3
    inputName = 'O'; % Firing rate
end

switch inputName
    case 'O'
        cirPos = exp(1i * NetPars.PrefStim/ NetPars.Width * pi);
        
        % BumpPos = squeeze(mean(bsxfun(@times, cirPos, O), 1)); % average over all neurons
        BumpPos = mean(bsxfun(@times, cirPos, Input), 1); % average over all neurons
        BumpPos = shiftdim(BumpPos, 1);
        BumpPos = BumpPos ./ abs(BumpPos); % normalize
        %         BumpPos(isnan(BumpPos)) = 0;
        IdxNaN = isnan(BumpPos);
%         BumpPos(IdxNaN(:)) = 360*rand(1, sum(IdxNaN(:)))-180;
        BumpPos(IdxNaN(:)) = exp(1i * (2*pi*rand(1, sum(IdxNaN(:))) - pi));
    case 'BumpPos'
        BumpPos = exp(1i * Input/ NetPars.Width * pi);
end

if nargout == 1
   BumpPos = angle(BumpPos)* NetPars.Width/pi; % angular value 
   return;
end

if isempty(BumpPos)
    % This case may happen when input variable is BumpPos
    szBumpPos = size(BumpPos);
    szBumpPos(2) = 1;
    meanBumpPos = nan(szBumpPos);
    varBumpPos = nan(szBumpPos);
    concBumpPos = nan(szBumpPos);
    mrlBumpPos = nan(szBumpPos);
else
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
    devBumpPos(devBumpPos>NetPars.Width) = devBumpPos(devBumpPos>NetPars.Width) - 2*NetPars.Width;
    varBumpPos = sum(devBumpPos.^2, 2) / (size(devBumpPos,2)-1);
end