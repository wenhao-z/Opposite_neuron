% Parameters for coupled two CANNs
% Wen-Hao Zhang, oct-06, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load default network parameters
defaultNetPars;

%% Specialized parameters of CANN

NetPars.numNets         = 2;
NetPars.numGroupPerNet  = 1;

% Switch into Gaussian connections; Default is vonMises
% NetPars.connFunc = 'Gaussian';
% NetPars.k         = 5e-4; % global inhibition strength
% NetPars.TunWidth  = 40; % Tuning width, the std. of tuning function. Unit: deg.

getDependentPars;

%% Caution: don't move below lines before parseNetPars
NetPars.Jrc         = (0.1: 0.1: 0.9)*NetPars.Jc;
NetPars.JrpRatio    = 0.1:0.1: 0.9; % same parameter for J12 and J21, relative to Jrc
NetPars.krpRatio    = 0.5; % Relative divisive normalization strength cross networks

% -----------------
% Input parameters
% -----------------
NetPars.Ampl        = (0.2: 0.1: 2)*NetPars.Uc;
NetPars.Ampl        = repmat(NetPars.Ampl, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

NetPars.Posi        = [NetPars.PrefStim(end/2-1), ...
    NetPars.PrefStim(end/2+1)]'; % the last element is used for different position

NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% Parameters of multiple trials
NetPars.tTrial     = 500 * NetPars.tau;
NetPars.nTrials    = 1; % number of trials

% NetPars.tLen = 50 * NetPars.tau;
% NetPars.nTrials = 20;

NetPars.tLen       = NetPars.nTrials * NetPars.tTrial;
NetPars.tStat      = 5 * NetPars.tau; % The starting time to make statistics

%%
% Parse network parameters
parseNetPars;
