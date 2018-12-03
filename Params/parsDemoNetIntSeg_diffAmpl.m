% Parameters for network of integration (congruent neuron) 
% and segregation (opposite neuron)
% Wen-Hao Zhang, Oct-09, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

parsDemoNetIntSeg;

% -----------------
% Input parameters
% -----------------
AmplRatio = 0.3: 0.1: 1.5;
[AmplRatio1, AmplRatio2] = meshgrid(AmplRatio, AmplRatio);
NetPars.AmplRatio = [AmplRatio1(:), AmplRatio2(:)]';
NetPars.AmplRatio   = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet, 1]);

NetPars.AmplBkg     = 1; 

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

% Parameters of multiple trials
NetPars.tLen    = 35 * NetPars.tau;
NetPars.nTrials = 25;
NetPars.tStat   = 15 * NetPars.tau; % The starting time to make statistics

%%
% Parse network parameters
parseNetPars;
