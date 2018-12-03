function savePath = getSavePath(Path_RootDir, NetPars)
% Generate the path for saving variables

% Wen-Hao Zhang, Oct-09, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

if (NetPars.numGroupPerNet == 1)
    switch NetPars.numNets
        case 1
            savePath = 'SingleCANN';
        case 2
            savePath = 'DoubleCANNs';
        otherwise
            savePath = sprintf('%dCoupledCANNs', NetPars.numNets);
    end
elseif (NetPars.numGroupPerNet == 2)
    savePath = 'NetIntSeg';
else
    savePath = [];
end

savePath = fullfile(Path_RootDir, 'Data', savePath);