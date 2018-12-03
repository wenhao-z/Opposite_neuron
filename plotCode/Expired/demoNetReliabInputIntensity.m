% Demonstration of network's reliability with input intensiy

% Wen-Hao Zhang, May-11, 2016


%% Net simulation
setPath;
% Load parameters
paramsNetDemo;

NetPars.Ampl = repmat(0.1: 0.1: 1, [2,1]) * NetPars.U0;
NetPars.tLen = 5 * NetPars.tau;
NetPars.NoiseFree = 0; % 1: no noise; 0: add noise;
NetPars.Posi = [0, 0]';
NetPars.cueCond = 0;

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

%% Produce random seeds
[seedIntNoisArray, seedExtNoisArray] = initRandSeed(NetPars, dimPar, size(parGrid));

%% Net Simulation
meanNetSim = cell(size(parGrid));
varNetSim = cell(size(parGrid));
concNetSim = cell(size(parGrid)); % concentration parameter of circular statistics

OHeight = cell(size(parGrid));
UHeight = cell(size(parGrid));
RateAvg = cell(size(parGrid));

% parpool(6);
tStart = clock;
for iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedIntNois = seedIntNoisArray(iterPar);
    netpars.seedExtNois = seedExtNoisArray(iterPar);
    
    if netpars.Jrc*(1 + netpars.JrpRatio) > netpars.Jc
        % Parameters of persistent activity
        BumpPos = nan(1,2);
        meanBumpPos = nan(4,1);
        varBumpPos = nan(4,1);
        concBumpPos = nan(4,1);
    else
        % Network input
        InputSet = makeNetInput(netpars.Posi, netpars);
        
        % Run simulation
        InputSet = simDecenNet2(InputSet, netpars);
        
        % Bump position
        O = InputSet.O(:,:, 2*NetPars.tau/ NetPars.dt+1 : end);
        [BumpPos, meanBumpPos, varBumpPos, concBumpPos] ...
            = statBumpPos(O, NetPars);
        
        RateAvg{iterPar} = mean(O,3);
        OHeight{iterPar} = shiftdim(max(mean(O, 3))', -ndims(parGrid));
        UHeight{iterPar} = shiftdim(max(mean(InputSet.U, 3))', -ndims(parGrid));
    end
    meanNetSim{iterPar} = shiftdim(meanBumpPos, -ndims(parGrid));
    varNetSim{iterPar} = shiftdim(varBumpPos, -ndims(parGrid));
    concNetSim{iterPar} = shiftdim(concBumpPos, -ndims(parGrid));
end

meanNetSim = cell2mat(meanNetSim); % the last dimension is index of network
varNetSim = cell2mat(varNetSim);
concNetSim = cell2mat(concNetSim);
OHeight = cell2mat(OHeight);
UHeight = cell2mat(UHeight);

tEnd = clock;

%%
figure
subplot(1,2,1);
plot([0, squeeze(OHeight(:,:,1))'], [0, squeeze(concNetSim(:,:,1))'] , '-o')
% hold on
% plot(squeeze(OHeight(:,:,1)), squeeze(concNetSim(:,:,3)) , '-or')
axis square;
xlabel('Peak firing rate (Hz)')
ylabel('Concentration')

subplot(1,2,2); 
plot(squeeze(OHeight(:,:,1)), squeeze(1./varNetSim(:,:,1)) , '-o')
axis square;