function [Rate, pctWin1] = simIntSegDMCircuit(NetPars, Iext)
% Simulate the network of two competing neurons (decision-making circuit)
% Wen-Hao Zhang, Mar-10, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Parameters for the decision-making circuit
% This part of code is in TEMPORAL!! REMOVE it in the final stage!

% NetPars.tau         = 1; % time constant
% NetPars.dt          = 0.01; % time step
% NetPars.Jint        = 0.7; % Mutual inibitory strength
% NetPars.nTrials     = 1;
% NetPars.kInt        = 1e-3;

%% Simulation
% Load fields from NetPars
tau = NetPars.tau;
dt = NetPars.dt;
% stdIntNois = NetPars.stdIntNois;
AmplBkg = NetPars.AmplBkg;
stdIntNois = sqrt(AmplBkg * NetPars.fanoFactor); % internal noise
kInt = NetPars.kInt;
% nTrials = NetPars.nTrials;
nTrials = size(Iext,3);

% Initialize matrices
Isyn = zeros(size(Iext));
Rate = zeros(size(Isyn));
IntNois = stdIntNois * randn(size(Iext));

Wint = NetPars.Jint*[0, -1; ...
    -1, 0]; % Mutual inhibition connection
pwr = 2;

for iterTrial = 1: nTrials
    for t = 1: length(Isyn)-1
        dI = -Isyn(:,t,iterTrial) + Iext(:,t,iterTrial) + Wint * Rate(:,t,iterTrial)...
            + IntNois(:,t,iterTrial)/sqrt(dt) + AmplBkg;
        Isyn(:,t+1,iterTrial) = Isyn(:,t,iterTrial) + dI *dt/tau;
        
        Irec = Isyn(:,t+1,iterTrial);
        Irec(Irec<0) = 0;
        Rate(:,t+1,iterTrial) = Irec.^pwr./ (1 + kInt*Irec.^pwr);
        
        % R = (Irec)./(1 - exp(-0.1*(Irec)));
        % R(isnan(R)) = 0;
        % Rate(:,t+1, iterTrial) = R;
    end
end

%% Make the statistics of winning probability
tLenStat = 10 * tau/dt-1;

flagWin1 = (squeeze(mean(Rate(1, end-tLenStat:end,:), 2)) > ...
    squeeze(mean(Rate(2, end-tLenStat:end,:), 2)));
pctWin1 = sum(flagWin1)/ nTrials; % The wining probability of 1st pool

%% Fold the output variables into a struct
% ......