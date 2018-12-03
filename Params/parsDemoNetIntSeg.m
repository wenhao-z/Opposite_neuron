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
NetPars.krpRatio    = 0.5; % Divsive normalization strength across groups, relative to NetPars.k

% getDependentPars;

%% Caution: don't move below lines before parseNetPars
NetPars.JrcRatio    = 0.3;  % Recurrent connection strength within the same network, relative to Jc
NetPars.JrpRatio    = 0.5; % same parameter for J12 and J21, relative to Jrc

% -----------------
% Input parameters
% -----------------
NetPars.AmplRatio   = 0.3:0.2:1.5; % relative to Uc
% NetPars.Ampl      = 1*NetPars.Uc;
NetPars.AmplRatio   = repmat(NetPars.AmplRatio, [NetPars.numNets * NetPars.numGroupPerNet, 1]);
NetPars.AmplBkg     = 1; 

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
% NetPars.Posi    = [0; 0];0.3
% NetPars.Posi    = [-30; 30];
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

NetPars.bAddNoise = 1;
NetPars.cueCond = 0: 2;

% NetPars.seedNois = sum(clock)*100;
NetPars.seedNois = 0;

% Parameters of multiple trials
NetPars.tLen    = 150 * NetPars.tau;
NetPars.nTrials = 1;
NetPars.tStat   = 10 * NetPars.tau; % The starting time to make statistics

%%
% Parse network parameters
parseNetPars;
