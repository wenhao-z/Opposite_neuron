function InputSet = simDecenNet_IntSeg(InputSet, NetPars, cueCond)
% A decentralized system for information integration and segregation
% The whole system is composed of several networks, with each is composed
% two groups of neurons: congruent and opposite neurons.
% Each neuronal group inside each network is modelled as a continuous attractor
% neural network

% (ref. W.H. Zhang. et al., JNS 2016 and W.H. Zhang et al., NIPS 2016)
% Each network module contains congruent (C) and opposite (O) neurons;
% Within each network, no recurrent connections between (C) and (O) neurons
% Across networks, (C) neurons are connected in a congruent manner; while
%                            (O) neurons are connected in an opposite manner.

% Author: Wen-Hao Zhang, Oct-6-2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

if nargin == 2
    cueCond = NetPars.cueCond;
end

% Unfold parameters from struct NetPars and InputSet
PrefStim    = NetPars.PrefStim;
Width       = NetPars.Width;
dt          = NetPars.dt;
tau         = NetPars.tau;
numNets     = NetPars.numNets;
Iext        = InputSet.Iext;
ExtNois     = InputSet.ExtNois;

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

%% Internal Nois
% set the seed of random number generator
if isfield(NetPars, 'seedIntNois')
    seedIntNois = NetPars.seedIntNois;
else
    seedIntNois = sum(clock*100);
end
s = RandStream('mt19937ar','Seed', seedIntNois);
RandStream.setGlobalStream(s);

if NetPars.bAddNoise
    switch cueCond
        case 1 % Cue 1
            Iext(:, 2:2:end, :) = zeros(size(Iext)./[1, 2,1]);
            ExtNois(:, 2:2:end, :) = zeros(size(Iext)./[1, 2,1]);
        case 2 % Cue 2
            Iext(:, 1:2:end, :) = zeros(size(Iext)./[1, 2,1]);
            ExtNois(:, 1:2:end, :) = zeros(size(Iext)./[1, 2,1]);
    end
    
    % the pseudo-noise sequence in 3 cueing conditions are the same
    IntNois = randn(size(Iext)) / sqrt(NetPars.dt);
    IntNois = bsxfun(@times, shiftdim(NetPars.stdIntNois, -1), IntNois);
    
    Iext = Iext + IntNois + ExtNois;
end

% Add the mean value of background inputs
Iext = Iext + NetPars.AmplBkg;

%% Iteration
% Initiation
U = zeros(size(Iext));
O = zeros(size(Iext));
pwr = 2; % the power of the u-r relation
divMat = repmat(diag(ones(1,numNets)), [1, numNets])';
nIterTrial = NetPars.tTrial/ NetPars.dt; % number of iteration for a trial

% O(:,1:2,1) = O(:,1:2,1) + 30;

for t = 1: size(Iext, 3) - 1
    OFt = O(:,:, t);
%     OFt(OFt < 10) = 0;
    OFt = fft(OFt);
    %     OFt = fft(O(:,:, t));
    
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
    if ~mod(t, nIterTrial)
        % Reset after finishing a trial
        U(:,:,t+1) = zeros(size(U(:,:,t+1)));
        O(:,:,t+1) = zeros(size(O(:,:,t+1)));
    else
        dU = (-U(:,:,t) + ISyn + Iext(:,:, t) ) * dt/tau;
        U(:,:,t+1) = U(:,:,t) + dU;
        
        % Synaptic input --> Firing rate
        Urec = U(:,:,t+1);
        Urec(Urec<0) = 0;
        Urec = Urec.^pwr;
        
        divU = sum(Urec, 1);
        divU = divU * divMat;
        divU = repmat(divU, [1, numNets]);
        
        O(:,:,t+1) = bsxfun(@rdivide, Urec, 1+ NetPars.k *divU);
    end
end

%% Fold variables into output struct
sizeArray = size(O);
sizeArray = [sizeArray(1:2), NetPars.tTrial/NetPars.dt, NetPars.nTrials];

InputSet.U = reshape(U, sizeArray);
InputSet.O = reshape(O, sizeArray);
InputSet.seedIntNois = seedIntNois;