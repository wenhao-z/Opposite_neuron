function InputSet = simDecenNet_Int(InputSet, NetPars, cueCond)
% A decentralized system consisting of several reciprocally connected networks.
% Each network is modelled as a continuous attractor neural network

% (ref. W.H. Zhang. et al., JNS 2016)
% Features:
% 1. Intra-network connections are translation-invariant bell shape function.
% 2. Inter-network connections have the same profile but different strength
%     with intra-network connections.
% 3. Each network has its own divisive normalization  pool

% INPUT
% InputSet.Iext: [N, numNets, Time]

% OUTPUT
% InputSet.U: [N, numNets, Time]; synaptic inputs
% InputSet.O: [N, nuvgxcmNets, Time]; firing rate

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
Iext        = InputSet.Iext;
ExtNois     = InputSet.ExtNois;

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
% make sure the kernel is strictly symmetric
% KerFt = (KerFt + flipud(KerFt))/2;
KerFt = fft(KerFt);

% Weight matrix [numNets, ,numNets]
% The matrix will be right multiply with neural firing rate R (N-by-2 array);
JMat = ones(NetPars.numNets) * NetPars.JrpRatio ...
    + diag(ones(1, NetPars.numNets)) * (1 - NetPars.JrpRatio);
JMat = JMat * NetPars.Jrc;

%% Internal Nois
% set the seed of random number generator
if isfield(NetPars, 'seedIntNois')
    seedIntNois = NetPars.seedIntNois;
else
    seedIntNois = sum(clock*100);
end
s = RandStream('mt19937ar','Seed', seedIntNois);
RandStream.setGlobalStream(s);

% This part needs to be revised to in accomodate with N>3 nets. (Oct-6, 2016)
if NetPars.bAddNoise
    switch cueCond
        case 1 % Cue 1
            Iext(:, 2, :) = zeros(size(Iext)./[1, 2,1]);
            ExtNois(:, 2, :) = zeros(size(Iext)./[1, 2,1]);
        case 2 % Cue 2
            Iext(:, 1, :) = zeros(size(Iext)./[1, 2,1]);
            ExtNois(:, 1, :) = zeros(size(Iext)./[1, 2,1]);
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
if ~isfield(NetPars, 'tTrial')
    NetPars.tTrial = NetPars.tLen;
    NetPars.nTrials = 1;
end
nIterTrial = NetPars.tTrial/ NetPars.dt; % number of iteration for a trial

for t = 1: size(Iext, 3) - 1
    OFt = fft(O(:,:, t));
    
    % Inputs received by congruent neurons
    ISyn = bsxfun(@times, KerFt,  OFt); % Nx2
    ISyn = ifft(ISyn) * JMat;
        
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
        
        O(:,:,t+1) = bsxfun(@rdivide, Urec, 1+ NetPars.k * divU);
    end
end

%% Fold variables into output struct
sizeArray = size(O);
sizeArray = [sizeArray(1:2), NetPars.tTrial/NetPars.dt, NetPars.nTrials];

InputSet.U = reshape(U, sizeArray);
InputSet.O = reshape(O, sizeArray);
InputSet.seedIntNois = seedIntNois;