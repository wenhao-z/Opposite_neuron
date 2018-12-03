% Default parameters of network
NetPars.N               = 180;  % The number of neurons
NetPars.numNets         = 2;  % number of networks
NetPars.numGroupPerNet  = 2; % number of group of neurons inside each network
%                              congruent and opposite neurons
NetPars.Width           = 180; % the range of parameter space from (-Width, Width), unit: deg

% Preferred stimulus of neurons (location on feature space)
PrefStim                = linspace(-NetPars.Width,NetPars.Width, NetPars.N+1)'; 
PrefStim(1)             = [];
NetPars.PrefStim        = PrefStim;
clear PrefStim

%% Temporal parameters
NetPars.tau  = 1; % Time constant of neuron activity
NetPars.tLen = 600 * NetPars.tau; % whole length of simulation
NetPars.dt   = NetPars.tau/100; % the iterative step

%% Connection
NetPars.connFunc = 'vonMises'; % or Gaussian
% NetPars.connFunc = 'Gaussian'; % or vonMises

switch NetPars.connFunc
    case 'Gaussian'
        NetPars.k         = 5e-4; % global inhibition strength
        NetPars.TunWidth  = 40; % Tuning width, the std. of tuning function. Unit: deg.
    case 'vonMises'
        NetPars.k           = 3e-4; % global inhibition strength
        NetPars.TunKappa    = 3; % Tuning width, concentration of von-Mises function, about 40 deg.
end
NetPars.krpRatio  = 0.5; % The inhibition strength from other groups of neurons within the same network. Relative to k
NetPars.k12Ratio  = 0;   % The inhibition strength from same type of neuron in different network. Relative to k
NetPars.JrcRatio = 0.5; % Recurrent connection strength within the same network, relative to Jc (the minimal recurrent 
%                         connection strength for the network to hold a persistent activity without feedforward inputs).
NetPars.JrpRatio = 0.1:0.1: 0.9; % Reciprocal connection strength between networks; same parameter for J12 and J21, relative to Jrc

%% Network input 
% -----------------------------
% Input intensity and location
% -----------------------------
% Peak intensity of feedforward inputs, [numNets*numGroupPerNet, 1]
NetPars.AmplRatio    = 0.2: 0.1: 1.6; % Relative to Uc, which is the persistent bump height without stimulus when Jrc = Jc 
NetPars.AmplRatio    = repmat(NetPars.AmplRatio, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% Intensity of background input
NetPars.AmplBkg = 1; 

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10: end)'; % every 20 degree
Posi            = [zeros(size(Posi)); Posi];
NetPars.Posi    = Posi;
NetPars.Posi    = repmat(NetPars.Posi, [NetPars.numGroupPerNet, 1]);
clear Posi

% ------------------
% Noise structure
% ------------------
NetPars.bAddNoise = 1; % 1: add noise; 0: noise free;
NetPars.PosiNoise = 0; % bool variable, 1: noise on the position of external input; 0: full noise

% Internal noise inside network
% The noise strength of all networks are the same for simplicity
NetPars.typeIntNois = 'Poisson';
switch NetPars.typeIntNois
    case 'Gaussian'
       NetPars.stdIntNois = [0: 10, 15: 5: 20]; % internal noise
    case 'Poisson'
        NetPars.fanoFactor = 0.5; % fano factor of noise
        NetPars.stdIntNois = sqrt(NetPars.AmplBkg * NetPars.fanoFactor); % internal noise
end
NetPars.stdIntNois = repmat(NetPars.stdIntNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% External noise associated with feedforward connections
NetPars.typeExtNois = 'Poisson'; % or 'Gaussian'
if strcmp(NetPars.typeExtNois, 'Gaussian')
    NetPars.stdExtNois = 0.3: 0.1: 0.6; % external noise
    NetPars.stdExtNois = repmat(NetPars.stdExtNois, [NetPars.numNets * NetPars.numGroupPerNet, 1]);
end

% ------------------
% Cueing conditions
% ------------------
NetPars.cueCond = 0: 2; % Cueing condition. 
%                         0: both cue; 
%                         1: only cue 1; 
%                         2: only cue 2.

% Random seed
NetPars.seedNois = 0;
% NetPars.seedIntNois = 0;
% NetPars.seedExtNois = sum(clock)*100;

NetPars.flagSeed = 1;
switch NetPars.flagSeed
    case 1
        NetPars.flagSeed = 'sameSeed';
        % use the same random seed for all parameters 
    case 2
        NetPars.flagSeed = 'SameSeedCueCond';
        % different random seed under different parameter settings, but for
        % each parameter set, the seeds under three cue conditions are
        % exactly the same
    case 3
        NetPars.flagSeed = 'diffSeed';
end

if NetPars.numGroupPerNet == 2
    NetPars.sameCORandSeed = 1; % 1: the feedforward inputs received by congruent
    %                                and opposite neurons have the same random seeds
    %                             0: different random seeds
end

NetPars = orderfields(NetPars);