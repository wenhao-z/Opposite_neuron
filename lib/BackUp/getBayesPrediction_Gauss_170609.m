function bayesRes = getBayesPrediction_Gauss(meanNetSim, varNetSim, dimPar)
% Calculate Bayesian predictions for integrated estimation
% Wen-hao Zhang, Apr-30-2015

%% Bayesian estimation of mean and variance

for iter = 1: length(dimPar)
    if strcmp(dimPar(iter).namePar, 'cueCond')
        dimCueCond = iter;
        break
    end
end
clear iter

expr = getExprCueCond('end-1:end', 0, ndims(varNetSim), dimCueCond);

meanNetSim1 = eval(['meanNetSim' expr]);
varNetSim1 = eval(['varNetSim' expr]);

varNetBayes = 1./  sum(1./varNetSim1, dimCueCond);

sz = ones(1, ndims(varNetBayes));
sz(dimCueCond) = 2;
coef = flip(varNetSim1, dimCueCond) ./...
    repmat(sum(varNetSim1, dimCueCond), sz);
meanNetBayes = sum(coef .* meanNetSim1, dimCueCond);
clear sz coef varNetSim1 meanNetSim1 expr

%% Actual  weights of direct cue for two nets
% Cue 1 is the direct cue for net 1, and cue 2 is the direct cue for net 2
expr01 = getExprCueCond(0, 1, ndims(varNetSim), dimCueCond); % combined cue condition, net 1
expr11 = getExprCueCond(1, 1, ndims(varNetSim), dimCueCond); % cue 1 condition, net 1
expr21 = getExprCueCond(2, 1, ndims(varNetSim), dimCueCond); % cue 2 condition, net 1

wDirCueNetSim1 = (eval(['meanNetSim', expr01]) - eval(['meanNetSim', expr21]) )./ ...
    (eval(['meanNetSim', expr11]) - eval(['meanNetSim', expr21]) );

expr02 = getExprCueCond(0, 2, ndims(varNetSim), dimCueCond); % combined cue condition, net 2
expr12 = getExprCueCond(1, 2, ndims(varNetSim), dimCueCond); % cue 1 condition, net 2
expr22 = getExprCueCond(2, 2, ndims(varNetSim), dimCueCond); % cue 2 condition, net 2

wDirCueNetSim2 = (eval(['meanNetSim', expr02]) - eval(['meanNetSim', expr12]) )./ ...
    (eval(['meanNetSim', expr22]) - eval(['meanNetSim', expr12]) );

wDirCueNetSim = cat(ndims(meanNetSim), wDirCueNetSim1, wDirCueNetSim2);
clear wDirCueNetSim1 wDirCueNetSim2

%% Predicted weights of direct cue for two nets
wDirCueNetBayes1 = eval(['varNetSim', expr21]) ./ (eval(['varNetSim', expr11]) + eval(['varNetSim', expr21]));
wDirCueNetBayes2 = eval(['varNetSim', expr12]) ./ (eval(['varNetSim', expr12]) + eval(['varNetSim', expr22]));
wDirCueNetBayes = cat(ndims(varNetSim), wDirCueNetBayes1, wDirCueNetBayes2);
clear wDirCueNetBayes1 wDirCueNetBayes2

%% Deviations of combination weight and variance
expr0 = getExprCueCond(0, 0, ndims(varNetSim), dimCueCond);
dMeanPerct = wDirCueNetSim - wDirCueNetBayes;
dVarPerct = eval(['varNetSim', expr0]) ./ varNetBayes - 1;

% Average deviations for direct cue of two coupled networks
dMeanPerctAvg = mean(dMeanPerct, ndims(meanNetSim));
dVarPerctAvg = mean(dVarPerct, ndims(varNetSim));


%% Fold into an output struct
bayesRes.varNetBayes = varNetBayes;
bayesRes.meanNetBayes = meanNetBayes;
bayesRes.dMeanPerct = dMeanPerct;
bayesRes.dVarPerct = dVarPerct;
bayesRes.dMeanPerctAvg = dMeanPerctAvg;
bayesRes.dVarPerctAvg = dVarPerctAvg;
bayesRes.wDirCueNetSim = wDirCueNetSim;
bayesRes.wDirCueNetBayes = wDirCueNetBayes;

%% embedded function 
    function expr = getExprCueCond(cueCond, netNum, nDIMs, dimCueCond)
        if isnumeric(cueCond)
            cueCond = num2str(cueCond+1);
        end
        expr = [repmat(':,', [1,dimCueCond-1]), cueCond ', ']; 
        
        % The last dim is network number. 0: for both nets
        if ~netNum
            expr = [expr, repmat(':,', [1, nDIMs-dimCueCond])];
            expr = expr(1:end-1);
        else
            expr = [expr, repmat(':,', [1, nDIMs-dimCueCond-1]), ...
                num2str(netNum)];
        end
        expr = ['(', expr, ')'];
    end

end