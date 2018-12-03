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

% Jc: the minimal recurrent connection strength for holding persistent
%     activities after switching off feedforward input
% Uc: the self-sustained bump height of synaptic input (without feedforward
%     inputs) when recurrent connection strength is Jc
% Oc: the self-sustained bump height of firing rate
switch NetPars.connFunc
    case 'Gaussian'
        NetPars.k         = 5e-4; % global inhibition strength
        NetPars.krpRatio  = 0.5; % The relative ratio of inhibition strength from other groups of neurons.
        NetPars.TunWidth  = 40; % Tuning width, the std. of tuning function. Unit: deg.
        NetPars.rho       = NetPars.N /(2*NetPars.Width); % The density of neurons
        
        NetPars.Jc  = sqrt(8*sqrt(2*pi)*(1+NetPars.krpRatio) * NetPars.k*NetPars.TunWidth/NetPars.rho);
        NetPars.Uc  = NetPars.Jc/(4*sqrt(pi)*NetPars.k*NetPars.TunWidth);
        NetPars.Oc  = NetPars.Uc^2 ...
            /(1 + sqrt(2*pi)*NetPars.k*NetPars.rho*NetPars.TunWidth*NetPars.Uc^2);  
    case 'vonMises'
        NetPars.k           = 3e-4; % global inhibition strength
        NetPars.krpRatio    = 0.5; % The relative ratio of inhibition strength from other groups of neurons.
        NetPars.TunKappa    = 3; % Tuning width, concentration of von-Mises function, about 40 deg.
        NetPars.rho         = NetPars.N /(2*pi); % The density of neurons per degree
        
        NetPars.Jc = 8*pi* (1+NetPars.krpRatio)* NetPars.k/ NetPars.rho;
        NetPars.Jc = NetPars.Jc * besseli(0, NetPars.TunKappa/2)^2 / besseli(0, NetPars.TunKappa);
        NetPars.Jc = sqrt(NetPars.Jc); % 0.93 is an empirical factor to cancel the approximation in theory!
        NetPars.Uc = NetPars.Jc*exp(NetPars.TunKappa/2)...
            / (4*pi*(1+NetPars.krpRatio)*NetPars.k * besseli(0, NetPars.TunKappa/2));
        NetPars.Oc = NetPars.Uc * exp(NetPars.TunKappa/2) * besseli(0, NetPars.TunKappa/2) ...
            /(NetPars.rho * NetPars.Jc * besseli(0, NetPars.TunKappa));       
end

NetPars.Jrc      = 0.5*NetPars.Jc; % Recurrent connection strength within the same network
NetPars.JrpRatio = 0.1:0.1: 0.9; % Reciprocal connection strength between networks; same parameter for J12 and J21, relative to Jrc

%% Network input 
% -----------------------------
% Input intensity and location
% -----------------------------
% Peak intensity of feedforward inputs, [numNets*numGroupPerNet, 1]
NetPars.Ampl    = (0.2: 0.1: 1.6)*NetPars.Uc; 
NetPars.Ampl    = repmat(NetPars.Ampl, [NetPars.numNets * NetPars.numGroupPerNet, 1]);

% Intensity of background input
NetPars.AmplBkg = 1; 

% Input location, [numNets*numGroupPerNet, 1]
% [1C, 2C, 1o, 2o]
Posi            = NetPars.PrefStim(end/2: 10:end)'; % every 20 degree
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