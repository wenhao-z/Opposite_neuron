function [InputSet, NetStat] = simDecenNet(InputSet, NetPars, outArgs)
% Simulate a decentralized neural network consisting of several
% reciprocally connected network modules.
% Each network module is modelled as a continuous attractor neural network
% (CANN)
% Each network may contain only congruent neurons, or both congruent and
% opposite neurons

% Author: Wen-Hao Zhang, June-4-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

if nargin == 2
    outArgs = struct('InputSet', [], 'NetStat', []);
end

switch NetPars.numGroupPerNet
    case 1
        % Each network only contains CONGRUENT neurons
        [InputSet, NetStat] = simDecenNet_Int(InputSet, NetPars, outArgs);
    case 2
        % Each network contains CONGRUENT and OPPOSITE neurons
        if isfield(NetPars, 'stdWrepNoisRatio') && (NetPars.stdWrepNoisRatio>0)
            [InputSet, NetStat] = simDecenNet_IntSeg_NoConv(InputSet, NetPars, outArgs);
        else
            [InputSet, NetStat] = simDecenNet_IntSeg(InputSet, NetPars, outArgs);
        end
end