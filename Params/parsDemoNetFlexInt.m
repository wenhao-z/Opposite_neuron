% Parameters for network of integration (congruent neuron) 
% and segregation (opposite neuron)
% Wen-Hao Zhang, Oct-09, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load default network parameters
defaultNetPars;

%% Specialized parameters of CANN

% Switch into Gaussian connections; Default is vonMises
% NetPars.connFunc = 'Gaussian';
% NetPars.k         = 5e-4; % global inhibition strength
% NetPars.TunWidth  = 40; % Tuning width, the std. of tuning function. Unit: deg.

% Set krpRatio BEFORE getDependentPars!!!
NetPars.kCORatio    = 0.5; % Divsive normalization strength across C/O groups within the same net.
NetPars.krpRatio    = 0.5; % Divisive normalization betweeen different networks

getDependentPars;

%% Caution: don't move below lines before parseNetPars
NetPars.Jrc         = 0.4*NetPars.Jc; % Recurrent connection weight within a network
NetPars.JrpRatio    = 0.5; % same parameter for J12 and J21, relative to Jrc
NetPars.bInt        = 1;   % Bool variable to indicate whether two nets cross-talk with each other
%                           1: integration; 0: segregation.

% -----------------
% Input parameters
% -----------------
NetPars.Ampl        = (0.3:0.2:1.5)*NetPars.Uc;  % Input intensity
% NetPars.Ampl        = 1*NetPars.Uc;
NetPars.Ampl        = repmat(NetPars.Ampl, [NetPars.numNets * NetPars.numGroupPerNet, 1]);
NetPars.AmplBkg     = 1; % Background input

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
% NetPars.Posi    = [0; 0];
% NetPars.Posi    = [-30; 30];
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

NetPars.bAddNoise = 1;
NetPars.cueCond   = 0: 2; % 0: cue1 + cue2; 1: only cue1; 2: only cue 2

NetPars.seedIntNois = rand * 1e4;
NetPars.seedExtNois = sum(clock)*100;

% Parameters of multiple trials
NetPars.tLen    = 150 * NetPars.tau;
NetPars.nTrials = 1;
NetPars.tStat   = 10 * NetPars.tau; % The starting time to make statistics

%% Parameters of Decision-making circuit
% Parameters for decision-making circuit
NetPars.Jint          = 0.8; % Mutual inibitory strength
NetPars.kInt          = 0.1;
NetPars.AmplBkgDM     = 30;
NetPars.fanoFactorDM  = 0.5;
NetPars.J_CO2DM       = 0.05; % raw experience: <= 0.1. (170507)
NetPars.J_biasInput   = 10; 

%%
% Parse network parameters
parseNetPars;
