% Plot the decoding resutls of single cue information from congruent and
% opposite neurons responses under combined cues, and also plot integration
% probablity.
% Wen-Hao Zhang, Nov-15, 2017
% wenhaoz1@andrew.cmu.edu

bSavePlot = 0;

setWorkPath;
datPath = fullfile(Path_RootDir, 'Data');
Folder = 'NetIntSeg';

% fileName = 'scanNetPars_170622_2134.mat'; % diff net pars.
fileName = 'scanNetPars_180611_2221.mat'; % sameSeedCueCond, multi-trial, new parameter range
fileName = 'scanNetPars_180613_2015.mat';
load(fullfile(datPath, Folder, fileName), ...
    'NetEstimRes', 'bayesRes', 'dimPar', 'NetPars');

%%
% Reference is the network's estimate under single cues
% Decoding from C and O neurons' response under combined conditions

namePar = {dimPar.namePar};
IdxCueCondDim = find(cellfun(@(x) strcmp(x, 'cueCond'), namePar));
IdxNeuronGroup = length(namePar) + 1;
IdxCueCombine = find(NetPars.cueCond==0);

% IdxPars.AmplRatio = 1:size(NetPars.AmplRatio,2);
IdxPars.AmplRatio = 1:8;
IdxPars.JrcRatio = 1:2;
IdxPars.JrpRatio = 1:7;
IdxPars.Posi = 1:10;
IdxPars.krpRatio = 1:3;
IdxPars.cueCond = 1:3; % combined condition

OHeightAvg = cropHighDimArray(NetEstimRes.OHeightAvg, namePar, IdxPars, [IdxNeuronGroup, IdxCueCondDim]);
meanBumpPos = cropHighDimArray(NetEstimRes.meanBumpPos, namePar, IdxPars, [IdxNeuronGroup, IdxCueCondDim]);
concBumpPos = cropHighDimArray(NetEstimRes.concBumpPos, namePar, IdxPars, [IdxNeuronGroup, IdxCueCondDim]);

% Permute bump position through adding random noise
bPermBumpPos = 1;
if bPermBumpPos
    
    szDat = size(meanBumpPos);
    randArray = repmat(rand(1,1,szDat(3)), szDat(1:2));
    
    randArray = (randArray - 0.5) * NetPars.Width*2;
    randArray = round(randArray);
    
    meanBumpPos = meanBumpPos + randArray;
    meanBumpPos = angle(exp(1i * meanBumpPos * pi/NetPars.Width))...
        * NetPars.Width/pi;
end

% Decoding the estimate under single-cue condition
vecNetRes = OHeightAvg .* exp(1i*meanBumpPos*pi/NetPars.Width);
% vecNetRes = concBumpPos .* exp(1i*meanBumpPos*pi/NetPars.Width);

decodeNetRes = cat(2, ...
    cat(1, ...
    (vecNetRes(1,IdxCueCombine,:) + vecNetRes(3,IdxCueCombine,:))/2, ... % p(s_1|x_1)
    (vecNetRes(2,IdxCueCombine,:) - vecNetRes(4,IdxCueCombine,:))/2), ...% p(s_2|x_1)
    cat(1, ...
    (vecNetRes(1,IdxCueCombine,:) - vecNetRes(3,IdxCueCombine,:))/2, ...% p(s_1|x_2)
    (vecNetRes(2,IdxCueCombine,:) + vecNetRes(4,IdxCueCombine,:))/2));% p(s_2|x_2)

meanDecodeNet = angle(decodeNetRes) * NetPars.Width/pi;
concDecodeNet = abs(decodeNetRes);

% Integration probability

% EstimResCrop = NetEstimRes;
% EstimResCrop.meanBumpPos = EstimResCrop.meanBumpPos(IdxPars.AmplRatio, IdxPars.JrcRatio, ...
%     IdxPars.JrpRatio, IdxPars.Posi, IdxPars.cueCond, IdxPars.krpRatio,:);
% EstimResCrop.concBumpPos = EstimResCrop.concBumpPos(IdxPars.AmplRatio, IdxPars.JrcRatio, ...
%     IdxPars.JrpRatio, IdxPars.Posi, IdxPars.cueCond, IdxPars.krpRatio,:);
% EstimResCrop.OHeightAvg = EstimResCrop.OHeightAvg(IdxPars.AmplRatio, IdxPars.JrcRatio, ...
%     IdxPars.JrpRatio, IdxPars.Posi, IdxPars.cueCond, IdxPars.krpRatio,:);
% EstimResCrop.mrlBumpPos = EstimResCrop.mrlBumpPos(IdxPars.AmplRatio, IdxPars.JrcRatio, ...
%     IdxPars.JrpRatio, IdxPars.Posi, IdxPars.cueCond, IdxPars.krpRatio,:);
% 
% [ProbInt1, ProbInt_Bayes] = estimIntProb(EstimResCrop, NetPars, dimPar, 1);
% ProbInt3 = estimIntProb(EstimResCrop, NetPars, dimPar, 3);

%% Plot decoding results
figure
hold on;
IdxCue1 = find(NetPars.cueCond == 1);
IdxCue2 = find(NetPars.cueCond == 2);

cSpec = lines(4);
sSpec = 'o';

% p(s_1|x_1)
hPlot(1) = plot(squeeze(meanBumpPos(1,IdxCue1,:)), squeeze(meanDecodeNet(1,IdxCue1-1,:)), ...
    sSpec, 'color', cSpec(1,:));
% p(s_1|x_2)
% hPlot(2) = plot(squeeze(meanBumpPos(1,IdxCue2,:)), squeeze(meanDecodeNet(1,IdxCue2-1,:)), ...
%     sSpec, 'color', cSpec(2,:));
plot([0,180], [0,180], '--k')
axis square;
axis([0 180 0 180])

% p(s_2|x_2)
hPlot(2) = plot(squeeze(meanBumpPos(2,IdxCue2,:)), squeeze(meanDecodeNet(2,IdxCue2-1,:)), ...
    sSpec, 'color', cSpec(3,:));
% p(s_2|x_1)
% plot(squeeze(meanBumpPos(2,IdxCue1,:)), squeeze(meanDecodeNet(2,IdxCue1-1,:)), ...
%     sSpec, 'color', cSpec(4,:))
plot(180*[-1,1], 180*[-1,1], '--k')

% legend('Mean under direct cue', 'Mean under indirect cue', 'location', 'northwest')
legend(hPlot, 's_1|x_1', 's_2|x_2')
xlabel('Actual mean of stimuli given single cue')
ylabel('Reconstructed mean of stimuli given single cues')
axis square;
axis([0 180 0 180])
set(gca, 'xtick', 0:45:180, 'ytick', 0:45:180)

% goodness of fit
X = [reshape(meanDecodeNet,1,[]), reshape(meanDecodeNet,1,[])]';
Y = [reshape(meanBumpPos(1:2,[IdxCue1, IdxCue2],:),1,[]), reshape(meanBumpPos(1:2,[IdxCue1, IdxCue2],:),1,[])]';

Idx = ((X(:) - Y(:)) > NetPars.Width);
X(Idx) = X(Idx) - 2*NetPars.Width;
Idx = ((X(:) - Y(:)) < -NetPars.Width);
X(Idx) = X(Idx) + 2*NetPars.Width;

[~, gof] = fit(X, Y, 'poly1');
axisLim = axis;
text(axisLim(2), axisLim(3), ['R^2=', num2str(gof.rsquare)], 'horizontalalignment', 'right', ...
    'verticalalignment', 'bottom');

title(...
{['AmplRatio=[', num2str(NetPars.AmplRatio(1,IdxPars.AmplRatio([1,end]))) , ...
'], JrcRatio=[' num2str(NetPars.JrcRatio(IdxPars.JrcRatio([1,end]))), ']' ], ...
['JrpRatio=[' num2str(NetPars.JrpRatio(IdxPars.JrpRatio([1,end]))), ...
'], krpRatio=[' num2str(NetPars.krpRatio(IdxPars.krpRatio([1,end]))), ']'  ]});

%% Plot integration probablity
subplot(1,2,1)
plot(ProbInt1(:), ProbInt_Bayes(:), '.')
axis square

subplot(1,2,2)
plot(ProbInt3(:), ProbInt_Bayes(:), '.')
axis square
