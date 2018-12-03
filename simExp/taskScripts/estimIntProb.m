function [ProbInt, ProbInt_Bayes, alpha] = estimIntProb(NetEstimRes, NetPars, dimPar, flagMethod)
% Read out integration probability from the joint activity of congruent and
% opposite neurons.
% Wen-Hao Zhang, Nov-15, 2017
% wenhaoz1@andrew.cmu.edu


dimCueCond = cellfun(@(x) strcmp(x, 'cueCond'), {dimPar.namePar});
dimCueCond = find(dimCueCond);

dimPosi = cellfun(@(x) strcmp(x, 'Posi'), {dimPar.namePar});
dimPosi = find(dimPosi);

switch flagMethod
    case 1
        % Using mean concentration to estimate integration probability
        concBumpPos = reorderField(NetEstimRes, 'concBumpPos');
        meanBumpPos = reorderField(NetEstimRes, 'meanBumpPos');
        szRes = size(concBumpPos);
        
        Numrator = concBumpPos(1,1,:).^2 - concBumpPos(1,3,:).^2; % k1 * k12 * cos(x_1-x_2)
        
        vecNetEstim = concBumpPos .* exp(1i * meanBumpPos*pi/NetPars.Width);
        kappa1 = abs((vecNetEstim(1,1,:) + vecNetEstim(1,3,:))/2);
        kappa2s = abs((vecNetEstim(1,1,:) - vecNetEstim(1,3,:))/2);
        kappaCos = Numrator ./ (kappa1 + kappa2s) / 4;
        kappaX = kappa1 .* kappa2s ./ (kappa1 + kappa2s);
        
        kappaCos = reshape(kappaCos, [szRes(3:end), 1]);
        kappaX = reshape(kappaX, [szRes(3:end), 1]);
    case 2
        % Using mean firing rate to estimate integration probability
        OHeight = NetEstimRes.OHeight; % [..., IdxGroup, time];
        BumpPos = NetEstimRes.BumpPos;
        
        % Get the data under combined cue condition,
        orderDim = [dimCueCond, ndims(OHeight)-1, ndims(OHeight), ...
            setdiff(1:ndims(OHeight)-2, dimCueCond)];
        OHeight = permute(OHeight, orderDim); % [cueCond, IdxGroup, time, Pars]
        BumpPos = permute(BumpPos, orderDim); % [cueCond, IdxGroup, time, Pars]
        szTmp = size(OHeight);
        OHeight = reshape(OHeight(1,:), szTmp(2:end)); % [IdxGroup, time, Pars]
        BumpPos = reshape(BumpPos(1,:), szTmp(2:end)); % [IdxGroup, time, Pars]
        
        % The size of parameters exlucding cueCond
        szNetPars = size(NetEstimRes.concBumpPos);
        szNetPars = szNetPars(1:end-1); % the last dim is IdxGroup
        szNetPars(dimCueCond) = [];
        
        % Estimate the concentration under two single-cue conditions.
        vecNetEstim = OHeight .* exp(1i * BumpPos*pi/NetPars.Width);
        decodeRes = cat(1, ...
            (vecNetEstim(1,:,:) + vecNetEstim(3,:,:))/2, ... % p(s_1|x_1)
            (vecNetEstim(1,:,:) - vecNetEstim(3,:,:))/2); % p(s_1|x_2), complex number, [cueCond, time, Pars]
        
        decodeRes = angle(decodeRes) * NetPars.Width/pi; % calculate the direction
        decodeRes = exp(1i * decodeRes/NetPars.Width*pi); % transform into complex representation
        
        mrlBumpPos = squeeze(abs(mean(decodeRes, 2))); % mean resultant length, average over time
        concBumpPos = mrl2Kappa(mrlBumpPos); % concentration parameter
        concBumpPos = reshape(concBumpPos,[2, szNetPars, 1]);
        
        concBumpPos_Comb = reorderField(NetEstimRes, 'concBumpPos');
        Numrator = concBumpPos_Comb(1,1,:).^2 - concBumpPos_Comb(1,3,:).^2; % k1 * k12 * cos(x_1-x_2)
        Numrator = reshape(Numrator, [1, szNetPars]);
        
        kappaCos = Numrator ./ (concBumpPos(1,:) + concBumpPos(2,:)) / 4;
        kappaX = 1./ sum(1./concBumpPos, 1);
    case 3
        % Use firing rate directly as concentration
        OHeightAvg = reorderField(NetEstimRes, 'OHeightAvg');
        meanBumpPos = reorderField(NetEstimRes, 'meanBumpPos');
        szRes = size(OHeightAvg);
        
        Numrator = OHeightAvg(1,1,:).^2 - OHeightAvg(1,3,:).^2; % k1 * k12 * cos(x_1-x_2)
        
        vecNetEstim = OHeightAvg .* exp(1i * meanBumpPos*pi/NetPars.Width);
        kappa1 = abs((vecNetEstim(1,1,:) + vecNetEstim(1,3,:))/2);
        kappa2s = abs((vecNetEstim(1,1,:) - vecNetEstim(1,3,:))/2);
        kappaCos = Numrator ./ (kappa1 + kappa2s) / 4;
        kappaX = kappa1 .* kappa2s ./ (kappa1 + kappa2s);
        
        kappaCos = reshape(kappaCos, [szRes(3:end), 1]);
        kappaX = reshape(kappaX, [szRes(3:end), 1]);
end

ProbInt = 1 - 1./(1 + exp(kappaCos - kappaX).*sqrt(2*pi * kappaX));

% Bayesian prediction of integration probability
if nargout >= 2
    mrlBumpPos = reorderField(NetEstimRes, 'mrlBumpPos');
    szRes = size(mrlBumpPos);
    
    kappaX_Bayes = mrlBumpPos(2,1,:) .* mrlBumpPos(3,1,:);
    kappaX_Bayes = squeeze(mrl2Kappa(kappaX_Bayes));
    kappaX_Bayes = reshape(kappaX_Bayes, [szRes(3:end), 1]);
    % kappaX_Bayes = (concBumpPos(2,1,:).^-1 + concBumpPos(3,1,:).^-1).^-1;
    
    dX = diff(NetPars.Posi(1:2,:),1,1);
    dX = shiftdim(dX(:), -dimPosi+1);
    
    ProbInt_Bayes = bsxfun(@times, kappaX_Bayes, cos(dX*pi/NetPars.Width)-1);
    ProbInt_Bayes = bsxfun(@times, exp(ProbInt_Bayes), sqrt(2*pi * kappaX_Bayes));
    
    % Find the ratio between the prior of two models to best explain
    % read out results
    % alpha: prior(int)/prior(seg)
    lsNetBayes = @(alpha) sum((1./(1+alpha * ProbInt_Bayes) - ProbInt).^2);
    
    alpha = fminbnd(lsNetBayes, 5, 30, ...
        optimset('TolX', 1e-15, 'MaxIter', 1e5));
    
    %     alpha = 5:0.01: 30;
    %     fval = lsNetBayes(alpha);
    %     [~, IdxMin] = min(fval);
    %     alpha = alpha(IdxMin);
    
    ProbInt_Bayes = 1 - 1./(1+ProbInt_Bayes * alpha);
end


    function varField = reorderField(NetEstimRes, varName)
        varField = NetEstimRes.(varName);
        % szEstim = size(varField);
        orderDim = [dimCueCond, ndims(varField), setdiff(1:ndims(varField)-1, dimCueCond)];
        
        varField = permute(varField, orderDim); % [cueCond, IdxGroup, ...]
        %         varField = reshape(varField, size(varField,1), size(varField,2), []); % [cueCond, IdxGroup, other parameters]
    end


end