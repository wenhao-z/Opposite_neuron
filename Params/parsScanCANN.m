% Parameters for a single CANN
% Wen-Hao Zhang, oct-06, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

%% Load default network parameters
defaultNetPars;

parseNetPars;

%% Specialized parameters of CANN
NetPars.numNets         = 1; 
NetPars.numGroupPerNet  = 1;

NetPars.Jrc     = (0.1: 0.1: 1.5)*NetPars.Jc;
NetPars.Ampl    = (0.1: 0.1: 2)*NetPars.Uc;
NetPars.Posi    = 0;

%% Parse network parameters
parseNetPars;