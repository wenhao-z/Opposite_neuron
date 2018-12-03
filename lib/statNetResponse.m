function NetStat = statNetResponse(BumpPos, OHeight, O, NetPars, NetStat)
% Get the statistics of network responses

% The input NetStat is a empty struct, the name of fields is used to
% indicate which variables to output.

% Author: Wen-Hao Zhang, June-7-2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

if isfield(NetStat, 'sumSqU')
    sumSqU = U;
    sumSqU(sumSqU<0) = 0;
    sumSqU = squeeze(sum(sumSqU.^pwr,1));
    NetStat.sumSqU = sumSqU;
end

if isfield(NetStat, 'OAvgXTime')
    NetStat.OAvgXTime = mean(O,3);
end

if isfield(NetStat, 'OStdXTime')
    NetStat.OStdXTime = std(O,[],3);
end

if isfield(NetStat, 'BumpPos')
   NetStat.BumpPos = BumpPos; 
end

if isfield(NetStat, 'OHeight')
   NetStat.OHeight = OHeight; 
end

if isfield(NetStat, 'OHeightAvg')
    OHeight = OHeight(:,NetPars.tStat/NetPars.dt+1:end, :);
    NetStat.OHeightAvg = mean(reshape(OHeight, size(OHeight, 1), []), 2);
    NetStat.OHeightAvg = squeeze(NetStat.OHeightAvg);
end

if isfield(NetStat, 'meanBumpPos')
    BumpPos = BumpPos(:,NetPars.tStat/NetPars.dt+1:end, :);
    BumpPos = reshape(BumpPos, size(BumpPos,1), []);
    [~, meanBumpPos, varBumpPos, concBumpPos, mrlBumpPos] = statBumpPos(BumpPos, NetPars, 'BumpPos');
    
    NetStat.meanBumpPos = meanBumpPos;
    NetStat.varBumpPos  = varBumpPos;
    NetStat.concBumpPos = concBumpPos;
    NetStat.mrlBumpPos  = mrlBumpPos; 
end