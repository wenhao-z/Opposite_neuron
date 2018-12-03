% Parameters for network of integration (congruent neuron) 
% and segregation (opposite neuron)
% Wen-Hao Zhang, Oct-09, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load default network parameters
defaultNetParsGate;

%% Specialized parameters of CANN

% Switch into Gaussian connections; Default is vonMises
% NetPars.connFunc = 'Gaussian';
% NetPars.k         = 5e-4; % global inhibition strength
% NetPars.TunWidth  = 40; % Tuning width, the std. of tuning function. Unit: deg.

getDependentPars;

%% Caution: don't move below lines before parseNetPars
NetPars.Jrc         = (0.3: 0.1: 0.6)*NetPars.Jc;
NetPars.JrpRatio    = 0.1:0.1: 0.9; % same parameter for J12 and J21, relative to Jrc

% NetPars.Jrc         = 0.5*NetPars.Jc;
% NetPars.JrpRatio    = 0.1:0.1:0.9; % same parameter for J12 and J21, relative to Jrc

% -----------------
% Input parameters
% -----------------
NetPars.Ampl        = (0.3: 0.2: 1.5)*NetPars.Uc;
% NetPars.Ampl        = (0.2: 0.1: 1.5)*NetPars.Uc;
% NetPars.Ampl        = 1*NetPars.Uc;
NetPars.Ampl        = repmat(NetPars.Ampl, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

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
% NetPars.tTrial     = 20 * NetPars.tau;
% NetPars.nTrials    = 25; % number of trials
% NetPars.tLen       = NetPars.nTrials * NetPars.tTrial;

NetPars.tLen = 500 * NetPars.tau;
NetPars.nTrials = 1;
NetPars.tStat = 50 * NetPars.tau;
%%
% Parse network parameters
parseNetPars;
