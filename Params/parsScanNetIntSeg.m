% Parameters for network of integration (congruent neuron) 
% and segregation (opposite neuron)
% Wen-Hao Zhang, Oct-09, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load default network parameters
defaultNetPars;

%% Specialized parameters of CANN
NetPars.krpRatio    = [0.05, 0.25, 0.5]; % Divsive normalization strength across groups within the same network. Relative to NetPars.k

%% Caution: don't move below lines before parseNetPars
NetPars.JrcRatio    = 0.3: 0.1: 0.4; % Recurrent connection strength within the same network, relative to Jc
% NetPars.JrcRatio    = 0.3: 0.1: 0.6; % Recurrent connection strength within the same network, relative to Jc
NetPars.JrpRatio    = 0.1: 0.1: 0.9; % same parameter for J12 and J21, relative to Jrc

% -----------------
% Input parameters
% -----------------
NetPars.AmplRatio   = 0.1: 0.2: 1.5; % relative to Uc
% NetPars.AmplRatio        = 0.2: 0.1: 1.5;
NetPars.AmplRatio   = repmat(NetPars.AmplRatio, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% Parameters of multiple trials
NetPars.tLen    = 35 * NetPars.tau; % 35
NetPars.nTrials = 25; % 25
NetPars.tStat   = 15 * NetPars.tau; % The starting time to make statistics 15

% Noises on connection matrix
NetPars.stdWrepNoisRatio = 0; %0.9; % Relative to the maximal weight of reciprocal connections
NetPars.seedWrepNois = sum(clock) * 100;

NetPars.flagSeed = 'SameSeedCueCond';

%% Parse network parameters
parseNetPars;
