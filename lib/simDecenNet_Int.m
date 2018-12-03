function [InputSet, NetStat] = simDecenNet_Int(InputSet, NetPars, outArgs)
% A decentralized system for information integration
% The whole system is composed of several networks, with each is modeled as
% a continuous attractor neural network (CANN)

% (ref. W.H. Zhang. et al., JNS 2016 and W.H. Zhang et al., NIPS 2016)
% Features:
% 1. Intra-network connections are translation-invariant bell shape function.
% 2. Inter-network connections have the same profile but different strength
%     with intra-network connections.

% Author: Wen-Hao Zhang, Mar-13-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Unfold parameters from struct NetPars and InputSet
PrefStim    = NetPars.PrefStim;
Width       = NetPars.Width;
dt          = NetPars.dt;
tau         = NetPars.tau;

if nargin == 2
    outArgs = struct('InputSet', [], 'NetStat', []);
end

%% Connection kernel with unit connection strength
switch NetPars.connFunc
    case 'Gaussian'
        TunWidth  = NetPars.TunWidth;
        KerFt = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
        KerFt = exp(-KerFt.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
    case 'vonMises'
        TunKappa  = NetPars.TunKappa;
        KerFt = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width)); % unit: rad
        KerFt = exp(TunKappa * cos(KerFt) )/(2*pi*besseli(0, TunKappa));
end
KerFt = fft(KerFt);

% Weight matrix
% The matrix will be right multiply with neural firing rate R (N-by-2 array);
JMat = [1, NetPars.JrpRatio; ...
    NetPars.JrpRatio, 1] * NetPars.Jrc;

% Divisive normalization
pwr = 2; % the power of the u-r relation
% Weight matrix for divisive normalization.
divMat = [1, NetPars.k12Ratio;...
    NetPars.k12Ratio, 1]*NetPars.k;

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
    end
    if NetPars.bAddNoise
        Iext = Iext + InputSet.IntNois + InputSet.ExtNois;
    end
    % Add the mean value of background inputs
    Iext = Iext + NetPars.AmplBkg;
    
    % -----------------------------------------------------------
    % Iteration over time
    for t = 1: size(Iext, 3) - 1
        OFt = fft(O(:,:, t));
        
        % Inputs received by congruent neurons
        ISyn = bsxfun(@times, KerFt,  OFt); % Nx2
        ISyn = ifft(ISyn) * JMat;
        
        % Update
        dU = (-U(:,:,t) + ISyn + Iext(:,:,t)) * dt/tau;
        U(:,:,t+1) = U(:,:,t) + dU;
        
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

