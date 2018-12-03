function NetPars = getDependentPars(NetPars)
% Calculate dependent parameters in NetPars.

% Wen-Hao Zhang, June-7, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Recalculate REFERENCE value for connection strength

% Jc: the minimal recurrent connection strength for holding persistent
%     activities after switching off feedforward input
% Uc: the self-sustained bump height of synaptic input (without feedforward
%     inputs) when recurrent connection strength is Jc
% Oc: the self-sustained bump height of firing rate

switch NetPars.connFunc
    case 'Gaussian'
        NetPars.rho = NetPars.N /(2*NetPars.Width); % The density of neurons
        NetPars.Jc  = sqrt(8*sqrt(2*pi)*(1+NetPars.krpRatio)* NetPars.k*NetPars.TunWidth/NetPars.rho);
        NetPars.Uc  = NetPars.Jc/(2*sqrt(pi)*(1+NetPars.krpRatio)*NetPars.k*NetPars.TunWidth);
        NetPars.Oc  = NetPars.Uc^2 ...
            /(1 + sqrt(2*pi)*NetPars.k*NetPars.rho*NetPars.TunWidth*NetPars.Uc^2);
    case 'vonMises'
        NetPars.rho = NetPars.N /(2*pi); % The density of neurons per degree
        NetPars.Jc  = 8*pi*(1+NetPars.krpRatio)* NetPars.k/ NetPars.rho;
        NetPars.Jc  = NetPars.Jc * besseli(0, NetPars.TunKappa/2)^2 / besseli(0, NetPars.TunKappa);
        NetPars.Jc  = sqrt(NetPars.Jc); % 0.93 is an empirical factor to cancel the approximation in theory!
        NetPars.Uc  = NetPars.Jc*exp(NetPars.TunKappa/2)...
            / (2*pi*(1+NetPars.krpRatio)*NetPars.k * besseli(0, NetPars.TunKappa/2));
        NetPars.Oc  = NetPars.Uc * exp(NetPars.TunKappa/2) * besseli(0, NetPars.TunKappa/2) ...
            /(NetPars.rho * NetPars.Jc * besseli(0, NetPars.TunKappa));
end

NetPars.Jrc = NetPars.JrcRatio * NetPars.Jc;
NetPars.Jrp = NetPars.JrpRatio * NetPars.Jrc;
NetPars.Ampl = NetPars.AmplRatio * NetPars.Uc;