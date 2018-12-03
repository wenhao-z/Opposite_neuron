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

getDependentPars;

%% Caution: don't move below lines before parseNetPars
NetPars.Jrc         = 0.4*NetPars.Jc;
NetPars.JrpRatio    = 0.5; % same parameter for J12 and J21, relative to Jrc

% -----------------
% Input parameters
% -----------------
NetPars.Ampl        = 3*NetPars.Uc;
NetPars.Ampl        = repmat(NetPars.Ampl, [NetPars.numNets * NetPars.numGroupPerNet, 1]);
NetPars.AmplBkg     = 1; 
NetPars.fanoFactor  = 0.5;

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
% NetPars.Posi    = [-30; 60];
% NetPars.Posi    = [-30; 30];
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% NetPars.tLen   = 600 * NetPars.tau; % whole length of simulation
% NetPars.nTrail = 50;
NetPars.bAddNoise = 1;
NetPars.cueCond = 0: 2;

NetPars.seedIntNois = rand * 1e4;
NetPars.seedExtNois = sum(clock)*100;

% Parameters of multiple trials
% NetPars.tTrial     = 50 * NetPars.tau;
% NetPars.nTrials    = 20; % number of trials
% NetPars.tLen       = NetPars.nTrials * NetPars.tTrial;

NetPars.tLen = 100 * NetPars.tau;
NetPars.nTrials = 1;
NetPars.tStat = 40 * NetPars.tau;
%%
% Parse network parameters
parseNetPars;
