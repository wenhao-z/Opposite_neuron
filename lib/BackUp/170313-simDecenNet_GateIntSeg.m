function [InputSet, NetStat] = simDecenNet_GateIntSeg(InputSet, NetPars, outArgs)
% A decentralized system for information integration and segregation
% The whole system is composed of several networks, with each is composed
% two groups of neurons: congruent (C) and opposite (O) neurons.
% C and O neurons share individual div. norm. pool.
% Each neuronal group inside each network is modelled as a continuous attractor
% neural network


% (ref. W.H. Zhang. et al., JNS 2016 and W.H. Zhang et al., NIPS 2016)
% Each network module contains congruent (C) and opposite (O) neurons;
% Within each network, no recurrent connections between (C) and (O) neurons
% Across networks, (C) neurons are connected in a congruent manner; 
%            while (O) neurons are connected in an opposite manner.

% Author: Wen-Hao Zhang, Mar-10-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Unfold parameters from struct NetPars and InputSet
PrefStim    = NetPars.PrefStim;
Width       = NetPars.Width;
dt          = NetPars.dt;
tau         = NetPars.tau;
% numNets     = NetPars.numNets;
% seedNois    = NetPars.seedNois;

if nargin == 2
    outArgs = struct('InputSet', [], 'NetStat', []);
end

%% Connection kernel with unit connection strength
switch NetPars.connFunc
    case 'Gaussian'
        TunWidth  = NetPars.TunWidth;
        KerFtCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
        KerFtOppo = angle(exp(1i*(PrefStim - PrefStim(1)+Width) *pi/Width))* Width/pi;
        KerFtCong = exp(-KerFtCong.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
        KerFtOppo = exp(-KerFtOppo.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
    case 'vonMises'
        TunKappa  = NetPars.TunKappa;
        KerFtCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width)); % unit: rad
        KerFtOppo = angle(exp(1i*(PrefStim - PrefStim(1)+Width) *pi/Width));
        KerFtCong = exp(TunKappa * cos(KerFtCong) )/(2*pi*besseli(0, TunKappa));
        KerFtOppo = exp(TunKappa * cos(KerFtOppo) )/(2*pi*besseli(0, TunKappa));
end
% make sure the kernel is strictly symmetric
% KerFtCong = (KerFtCong + flipud(KerFtCong))/2;
% KerFtOppo = (KerFtOppo + flipud(KerFtOppo))/2;

KerFtCong = fft(KerFtCong);
KerFtOppo = fft(KerFtOppo);

% Weight matrix
% The matrix will be right multiply with neural firing rate R (N-by-2 array);
JMat = [1, NetPars.JrpRatio; ...
    NetPars.JrpRatio, 1] * NetPars.Jrc;
JMatOppo = [diag(diag(JMat)); JMat - diag(diag(JMat))];

%% Initiation
% [N, numNets, Time, nTrials]
sizeU = [NetPars.N, NetPars.numNets*NetPars.numGroupPerNet, ...
    NetPars.tLen/NetPars.dt, NetPars.nTrials];
if isfield(outArgs, 'InputSet')
    UArray = zeros(sizeU);
    OArray = zeros(sizeU);
end
if isfield(outArgs, 'NetStat')
    if isfield(outArgs.NetStat, 'BumpPos')
        BumpPos = zeros(sizeU(2:end));
    end
    if isfield(outArgs.NetStat, 'OHeight')
        OHeight = zeros(size(BumpPos));
    end
end

pwr = 2; % the power of the u-r relation
% divMat = repmat(diag(ones(1,numNets)), [1, numNets])';

%% Iteration
for iterTrial = 1: NetPars.nTrials
    U = zeros(sizeU(1:3));
    O = zeros(sizeU(1:3));
    
    % Generate new noise sequence of every trial
    InputSet = makeNetInput(InputSet, NetPars, ...
        struct('ExtNois', [], 'IntNois', []));
    
    %     [InputSet, NetPars] = makeNetNoise(InputSet, NetPars);
    if NetPars.bAddNoise
        Iext = InputSet.Iext + InputSet.IntNois + InputSet.ExtNois;
    end
    % Add the mean value of background inputs
    Iext = Iext + NetPars.AmplBkg;
    
    % -----------------------------------------------------------
    % Iteration over time
    for t = 1: size(Iext, 3) - 1
        OFt = fft(O(:,:, t));
        
        % Inputs received by congruent neurons
        ICong = bsxfun(@times, KerFtCong,  OFt(:, 1:end/2)); % Nx2
        ICong = ifft(ICong) * JMat;
        
        % Inputs received by opposite neurons
        IOppo = [bsxfun(@times, KerFtCong, OFt(:,end/2+1:end) ), ...
            bsxfun(@times, KerFtOppo, OFt(:,end/2+1:end)) ];
        % IOppo(:,1:2): recurrent inputs (congruent connection)
        % IOppo(:,3:4): reciprocal inputs (opposite connection)
        IOppo = ifft(IOppo* JMatOppo);
        ISyn = [ICong, IOppo];
        
        % Update
        dU = (-U(:,:,t) + ISyn + Iext(:,:,t) ) * dt/tau;
        U(:,:,t+1) = U(:,:,t) + dU;
        
        % Synaptic input --> Firing rate
        Urec = U(:,:,t+1);
        Urec(Urec<0) = 0;
        divU = sum(Urec.^pwr, 1);
        %         divU = divU * divMat;
        %         divU = repmat(divU, [1, numNets]);
        
        O(:,:,t+1) = bsxfun(@rdivide, Urec, 1+ NetPars.k *divU);
        
        % ----------------------------------------------
        % Dynamics for the Int and Seg interneurons, Mar-10, 2017
        % Let me first test the network outside(4:04pm, Mar-10,2017)
        %         Oint
        %         Oseg
        
        
    end
    
    if isfield(outArgs, 'InputSet')
        UArray(:,:,:, iterTrial) = U;
        OArray(:,:,:, iterTrial) = O;
    end
    
    % Make statistics of network's activities
    % Calculate the bump position and height
    if exist('BumpPos', 'var')
        BumpPos(:,:,iterTrial) = statBumpPos(O, NetPars);
    end
    if exist('OHeight', 'var')
        OHeight(:,:,iterTrial) = sum(O, 1)*exp(NetPars.TunKappa) ...
            /(NetPars.rho*2*pi*besseli(0, NetPars.TunKappa)); % Average over neurons
    end
end

%% Estimate the statistics of network activities
% Q: how do I judge whether winner-take-all happens? Dec-31, 2016

if isfield(outArgs, 'NetStat')
    % Get the winner of groups
    tLenStat = 10 * NetPars.tau/NetPars.dt-1;
    flagCWin = (squeeze(mean(OHeight(1, end-tLenStat:end, :),2)) > ...
        squeeze(mean(OHeight(3, end-tLenStat:end, :),2)));
    pctCWin = mean(flagCWin);
    pctCWin = kron([pctCWin; 1-pctCWin], ones(2,1));
    
    % Calculate the statistics according to winner of groups
    BumpPosTmp = reshape(BumpPos(:,NetPars.tStat/NetPars.dt+1:end, flagCWin), sizeU(2), []);
    [~, meanBumpPos, varBumpPos, concBumpPos, mrlBumpPos] = statBumpPos(BumpPosTmp, NetPars, 'BumpPos');
    BumpPosTmp = reshape(BumpPos(:,NetPars.tStat/NetPars.dt+1:end, ~flagCWin), sizeU(2), []);
    [~, meanBumpPosOp, varBumpPosOp, concBumpPosOp, mrlBumpPosOp] = statBumpPos(BumpPosTmp, NetPars, 'BumpPos');
    clear BumpPosTmp
    
    OHeightAvg = mean(reshape(OHeight(:,NetPars.tStat/NetPars.dt+1:end, flagCWin), sizeU(2), []), 2);
    OHeightAvg = squeeze(OHeightAvg);
    
    % The output fields are dependent on the struct defined outside the
    % function. Dec-7, 2017
    for varName = fieldnames(outArgs.NetStat)'
        NetStat.(varName{1}) = eval(varName{1});
    end
end

%% Fold variables into output struct
if isfield(outArgs, 'InputSet')
    InputSet.U = UArray;
    InputSet.O = OArray;
%     InputSet = struct('U', UArray, 'O', OArray);
end

