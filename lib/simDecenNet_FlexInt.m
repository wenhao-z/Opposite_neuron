function [InputSet, NetStat] = simDecenNet_FlexInt(InputSet, NetPars, outArgs)
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

% Author: Wen-Hao Zhang, May-01-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Unfold parameters from struct NetPars and InputSet
PrefStim    = NetPars.PrefStim;
Width       = NetPars.Width;
dt          = NetPars.dt;
tau         = NetPars.tau;
numNets     = NetPars.numNets;
rho         = NetPars.rho;
% seedNois    = NetPars.seedNois;

% Parameters of decision making circuit
AmplBkgDM   = NetPars.AmplBkgDM;

if nargin == 2
    outArgs = struct('InputSet', [], 'NetStat', []);
end

if ~NetPars.bInt
    NetPars.JrpRatio = 0;
    NetPars.krpRatio = 0;    
end
%% Connection kernel with unit connection strength
dPrefStim = mean(diff(PrefStim));
switch NetPars.connFunc
    case 'Gaussian'
        TunWidth  = NetPars.TunWidth;
        KerFtCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
        KerFtOppo = angle(exp(1i*(PrefStim - PrefStim(1)+dPrefStim+Width) *pi/Width))* Width/pi;
        KerFtCong = exp(-KerFtCong.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
        KerFtOppo = exp(-KerFtOppo.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);
    case 'vonMises'
        TunKappa  = NetPars.TunKappa;
        KerFtCong = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width)); % unit: rad
        KerFtOppo = angle(exp(1i*(PrefStim - PrefStim(1)+dPrefStim+Width) *pi/Width));
        KerFtCong = exp(TunKappa * cos(KerFtCong) )/(2*pi*besseli(0, TunKappa));
        KerFtOppo = exp(TunKappa * cos(KerFtOppo) )/(2*pi*besseli(0, TunKappa));
end
KerFtCong = fft(KerFtCong);
KerFtOppo = fft(KerFtOppo);

% Weight matrix
% The matrix will be right multiply with neural firing rate R (N-by-2 array);
JMat = [1, NetPars.JrpRatio; ...
    NetPars.JrpRatio, 1] * NetPars.Jrc;
JMatOppo = [diag(diag(JMat)); JMat - diag(diag(JMat))];

% ---------------------------
% Divisive normalization
pwr = 2; % the power of the u-r relation
% Weight matrix for divisive normalization.
divMat = [1, NetPars.kCORatio; ...
    NetPars.kCORatio, 1];
divMat = kron(divMat, eye(numNets));

% ---------------------------------------------------------------
% Weight matrix for mutual inhibition in decision making circuit
Wint = NetPars.Jint*[0, -1; ...
    -1, 0]; % Mutual inhibition connection

%% Initiation
% [N, numNets, Time, nTrials]
sizeU = [NetPars.N, NetPars.numNets*NetPars.numGroupPerNet, ...
    NetPars.tLen/NetPars.dt, NetPars.nTrials];
if isfield(outArgs, 'InputSet')
    UArray = zeros(sizeU);
    OArray = zeros(sizeU);
    rDMArray = zeros([2, sizeU(3:4)]);
    InputDMArray = zeros([2,sizeU(3:4)]);
end
if isfield(outArgs, 'NetStat')
    BumpPos = zeros(sizeU(2:end));
    OHeight = zeros(size(BumpPos));
end

%% Iteration
for iterTrial = 1: NetPars.nTrials
    U = zeros(sizeU(1:3));
    O = zeros(sizeU(1:3));
    InputDM = zeros([2,sizeU(3)]);
    divUArray = zeros([4,sizeU(3)]);
    uDM = zeros([2,sizeU(3)]);
    rDM = zeros([2,sizeU(3)]);
    uDM(:,1) = NetPars.AmplBkgDM;
    rDM(:,1) = uDM(:,1).^pwr./ (1 + NetPars.kInt*uDM(:,1).^pwr);
    biasInput = 0;
    
    % ------------------------------------------
    % Generate new noise sequence of every trial
    InputSet = makeNetInput(InputSet, NetPars, ...
        struct('ExtNois', [], 'IntNois', []));
    
    Iext = InputSet.Iext;
    if NetPars.bAddNoise
        Iext = Iext + InputSet.IntNois + InputSet.ExtNois;
    end
    % Add the mean value of background inputs
    Iext = Iext + NetPars.AmplBkg;
    
    % Background noise fed to decision-making circuit
    IntNoisDM = sqrt(NetPars.AmplBkg * NetPars.fanoFactorDM) ...
        * randn([2, sizeU(3)])/sqrt(dt);
        
    % -----------------------------------------------------------
    % Iteration over time
    for t = 1: size(Iext, 3) - 1
        OFt = fft(O(:,:, t));
        
        % Inputs received by congruent neurons
        ICong = bsxfun(@times, KerFtCong,  OFt(:, 1:end/2)); % Nx2
        ICong = ifft(ICong);
        IrpCongHeight = sum(ICong * JMat(2), 1)*exp(NetPars.TunKappa/2) ...
            /(NetPars.rho*2*pi*besseli(0, NetPars.TunKappa/2)); % Average over neurons
        
        % Inputs received by opposite neurons
        IOppo = [bsxfun(@times, KerFtCong, OFt(:,end/2+1:end) ), ...
            bsxfun(@times, KerFtOppo, OFt(:,end/2+1:end))];
        % IOppo(:,1:2): recurrent inputs (congruent connection)
        % IOppo(:,3:4): reciprocal inputs (opposite connection)
        IOppo = ifft(IOppo);
        
        divIrp = [ICong*JMat(2), IOppo(:,3:4)*JMatOppo(2)];
        
        % Total synaptic input
        ISyn = [ICong*JMat, IOppo*JMatOppo];
        
        % -----------------------------------
        % Update
        dU = (-U(:,:,t) + ISyn + Iext(:,:,t)) * dt/tau;
        U(:,:,t+1) = U(:,:,t) + dU;
        
        % -----------------------------------
        % Synaptic input --> Firing rate
        Urec = U(:,:,t+1);
        Urec(Urec<0) = 0;
        Urec = Urec.^pwr;
        divU = sum(Urec, 1);
        
        if NetPars.bInt
            % Inter-module inputs targeting on div. norm neuron
            divIrp(divIrp<0) = 0;
            divIrp = sum(divIrp.^pwr, 1);
            divU = divU + divIrp * NetPars.krpRatio;
        end
        InputDM(:,t) = [mean(divU(1:2)); mean(divU(3:4))]; % Feedforward input to DM circuit       
        divU = divU * divMat;
        divUArray(:,t) = divU';
        
        O(:,:,t+1) = bsxfun(@rdivide, Urec, 1+ NetPars.k*divU);
                
        % ---------------------------------------------------------------
        % Decision making circuit
        % Feedforward input from div. norm neurons affiliated with C/O
        % neurons to decision making neurons
        kappaX = mean(IrpCongHeight) + NetPars.Ampl(1);
        % biasInput = kappaX - log(2*pi*kappaX)/2;
        % biasInput(isinf(biasInput)) = 0;
        biasInput = biasInput + (-biasInput + kappaX - log(2*pi*kappaX)/2) * dt/tau;
        
        InputDM(:,t) = InputDM(:,t)/ (rho * exp(-TunKappa) * 2*pi*besseli(0, TunKappa));
        InputDM(2,t) = InputDM(2,t) + NetPars.J_biasInput*biasInput';
        % InputDM(1,t) = InputDM(1,t) - NetPars.J_biasInput*biasInput';
        InputDM(:,t) = -flipud(InputDM(:,t))*NetPars.J_CO2DM;
        
        % Iteration of decision-making circuit
        duDM = -uDM(:,t) + InputDM(:,t) + Wint * rDM(:,t)...
            + IntNoisDM(:,t) + AmplBkgDM;
        uDM(:,t+1) = uDM(:,t) + duDM*dt/tau;
        
        uDMRec = uDM(:,t+1);
        uDMRec(uDMRec<0) = 0;
        rDM(:,t+1) = uDMRec.^pwr./ (1 + NetPars.kInt*uDMRec.^pwr);
    end
    
    if isfield(outArgs, 'InputSet')
        UArray(:,:,:, iterTrial) = U;
        OArray(:,:,:, iterTrial) = O;
        rDMArray(:,:,iterTrial) = rDM;
        InputDMArray(:,:,iterTrial) = InputDM;
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

if isfield(outArgs.NetStat, 'OAvgXTime')
    OAvgXTime = mean(O,3);
end

if isfield(outArgs, 'NetStat')
    % Estimate the mean and concentration of bump position
    BumpPosTmp = reshape(BumpPos(:,NetPars.tStat/NetPars.dt+1:end, :), sizeU(2), []);
    [~, meanBumpPos, varBumpPos, concBumpPos, mrlBumpPos] = statBumpPos(BumpPosTmp, NetPars, 'BumpPos');
    
    % Time averaged bump height
    OHeightAvg = mean(reshape(OHeight(:,NetPars.tStat/NetPars.dt+1:end, :), sizeU(2), []), 2);
    OHeightAvg = squeeze(OHeightAvg);
    
    % Time averaged activity of divisive normalization neuron
    divUAvg = mean(divUArray, 2);
    
    % Time averaged input received by decision-making neurons
    InputDMAvg = mean(InputDM, 2);
    
    % Winning probability of decision making circuit
    tLenStat = 10 * tau/dt-1;
    flagWin1 = (squeeze(mean(rDMArray(1, end-tLenStat:end,:), 2)) > ...
        squeeze(mean(rDMArray(2, end-tLenStat:end,:), 2)));
    pctWin1 = sum(flagWin1)/ NetPars.nTrials; % The wining probability of 1st pool

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
    InputSet.rDM = rDMArray;
    InputSet.InputDM = InputDMArray;
    InputSet.divU = divUArray;
end

