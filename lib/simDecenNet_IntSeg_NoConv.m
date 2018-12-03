function [InputSet, NetStat] = simDecenNet_IntSeg_NoConv(InputSet, NetPars, outArgs)
% A decentralized system for information integration and segregation
% The whole system is composed of several networks, with each is composed
% two groups of neurons: congruent and opposite neurons.
% Each neuronal group inside each network is modelled as a continuous attractor
% neural network

% (ref. W.H. Zhang. et al., JNS 2016 and W.H. Zhang et al., NIPS 2016)
% Each network module contains congruent (C) and opposite (O) neurons;
% Within each network, no recurrent connections between (C) and (O) neurons
% Across networks, (C) neurons are connected in a congruent manner; while
%                  (O) neurons are connected in an opposite manner.

% Author: Wen-Hao Zhang, July-05-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Unfold parameters from struct NetPars and InputSet
PrefStim    = NetPars.PrefStim;
Width       = NetPars.Width;
dt          = NetPars.dt;
tau         = NetPars.tau;
numNets     = NetPars.numNets;

if nargin == 2
    outArgs = struct('InputSet', [], 'NetStat', []);
end

%% Connection matrix
% The matrix will be right multiply with neural firing rate R (N-by-2 array);
JMat = [1, NetPars.JrpRatio; ...
    NetPars.JrpRatio, 1] * NetPars.Jrc;

dPrefStim = mean(diff(PrefStim));
switch NetPars.connFunc
    case 'Gaussian'
        TunWidth  = NetPars.TunWidth;
        WCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
        WOppo = angle(exp(1i*(PrefStim - PrefStim(1)+dPrefStim+Width) *pi/Width))* Width/pi;
        WCong = exp(-WCong.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
        WOppo = exp(-WOppo.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
    case 'vonMises'
        TunKappa  = NetPars.TunKappa;
        WCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width)); % unit: rad
        WOppo = angle(exp(1i*(PrefStim - PrefStim(1)+dPrefStim+Width) *pi/Width));
        WCong = exp(TunKappa * cos(WCong) )/(2*pi*besseli(0, TunKappa));
        WOppo = exp(TunKappa * cos(WOppo) )/(2*pi*besseli(0, TunKappa));
end
WCong = gallery('circul',WCong);
WOppo = gallery('circul',WOppo);

WOppo = [JMat(1) * WCong, JMat(3) * WOppo; ...
    JMat(2) * WOppo, JMat(4) * WCong];
WCong = kron(JMat, WCong);
% WOppo = kron(JMat, WOppo);

% ---------------------------------------------------------------
% Add noises onto reciprocal connections
WrpMax = max(WCong(1,:));

rng(NetPars.seedWrepNois);
Nois_WCong = NetPars.stdWrepNoisRatio * WrpMax * randn(size(WCong));
% Nois_WCong = NetPars.stdWrepNoisRatio .* sqrt(WCong) .* randn(size(WCong));

% Nois_WCong(1:end/2,1:end/2) = 0;
% Nois_WCong(end/2+1:end,end/2+1:end) = 0;

% Nois_WOppo = NetPars.stdWrepNoisRatio * WrpMax * randn(size(WOppo));
% Nois_WOppo = NetPars.stdWrepNoisRatio .* sqrt(WOppo) .* randn(size(WOppo));
% Nois_WOppo(1:end/2,1:end/2) = 0;
% Nois_WOppo(end/2+1:end,end/2+1:end) = 0;

Nois_WOppo = Nois_WCong;
Nois_WOppo(1:end/2, end/2+1:end) = Nois_WOppo(1:end/2, [3*end/4+1:end, end/2+1:3*end/4]);
Nois_WOppo(end/2+1:end, 1:end/2) = Nois_WOppo(end/2+1:end, [end/4+1:end/2, 1:end/4]);

WCong = WCong + Nois_WCong;
WOppo = WOppo + Nois_WOppo;

WCong(WCong<0) = 0;
WOppo(WOppo<0) = 0;

% ---------------------------------------------------------------
% Divisive normalization
pwr = 2; % the power of the u-r relation
% Weight matrix for divisive normalization.
divMat = [1, NetPars.krpRatio; ...
    NetPars.krpRatio, 1] * NetPars.k;
divMat = kron(divMat, eye(numNets));

%% Initiation
% [N, numNets, Time, nTrials]
sizeU = [NetPars.N, NetPars.numNets*NetPars.numGroupPerNet, ...
    NetPars.tLen/NetPars.dt, NetPars.nTrials];
if isfield(outArgs, 'InputSet')
    UArray = zeros(sizeU);
    OArray = zeros(sizeU);
end
if isfield(outArgs, 'NetStat')
    BumpPos = zeros(sizeU(2:end));
    OHeight = zeros(size(BumpPos));
end

%% Iteration
for iterTrial = 1: NetPars.nTrials
    U = zeros(sizeU(1:3));
    O = zeros(sizeU(1:3));
    
    % ------------------------------------------
    % Generate new noise sequence of every trial
    InputSet = makeNetInput(InputSet, NetPars, ...
        struct('ExtNois', [], 'IntNois', []));
    
    if size(InputSet.Iext, 3) == 1
        Iext = repmat(InputSet.Iext, [1,1,InputSet.szIext(3)]);
    else
        Iext = InputSet.Iext;
    end
    if NetPars.bAddNoise
        Iext = Iext + InputSet.IntNois + InputSet.ExtNois;
    end
    % Add the mean value of background inputs
    Iext = Iext + NetPars.AmplBkg;
        
    % -----------------------------------------------------------
    % Iteration over time
    for t = 1: size(Iext, 3) - 1
        % Inputs received by congruent neurons
        ICong = WCong * reshape(O(:,1:end/2, t), [],1);
        
        % Inputs received by opposite neurons
        IOppo = WOppo * reshape(O(:,end/2+1:end, t), [],1);
        
        % Total synaptic input
        ISyn = [reshape(ICong, [], 2), reshape(IOppo, [], 2)];
        
        % ---------------------------------------------------------
        % Update
        dU = (-U(:,:,t) + ISyn + Iext(:,:,t)) * dt/tau;
        U(:,:,t+1) = U(:,:,t) + dU;
        
        % -------------------------------
        % Synaptic input --> Firing rate
        Urec = U(:,:,t+1);
        Urec(Urec<0) = 0;
        Urec = Urec.^pwr;
        divU = sum(Urec, 1);
        divU = divU * divMat;
        
        O(:,:,t+1) = bsxfun(@rdivide, Urec, 1+divU);
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
if isfield(outArgs, 'NetStat')
    NetStat = statNetResponse(BumpPos, OHeight, O, NetPars, outArgs.NetStat);
end

%% Fold variables into output struct
if isfield(outArgs, 'InputSet')
    InputSet.U = UArray;
    InputSet.O = OArray;
end

