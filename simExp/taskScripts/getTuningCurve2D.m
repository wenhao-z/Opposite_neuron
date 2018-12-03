% Two-dim tuning curves of neurons in a decentralized network

% Author: Wen-Hao Zhang, June-2, 2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Parameters of inputs
NetPars.cueCond = [1,2,0]; % cue 1, cue 2, both cues
NetPars.tLen = 20 * NetPars.tau;

Posi = -180: 45: 180;
PosiArray = zeros([length(Posi)*ones(1,2),2]);
[PosiArray(:,:,1), PosiArray(:,:,2)] = meshgrid(Posi, Posi);
PosiArray = permute(PosiArray, [3,1,2]);
NetPars.Posi = reshape(PosiArray, 2, []);
NetPars.Posi = repmat(NetPars.Posi, [NetPars.numGroupPerNet,1,1]);
NetPars.flagSeed = 'sameSeed';

% NetPars.AmplRatio = 1.2*ones(2, 3);
% NetPars.AmplRatio(1,:) = [0.7, 0.35, 0.15];
NetPars.AmplRatio = 1.5*ones(2, 2);
NetPars.AmplRatio(1,:) = [0.35, 0.1];
NetPars.AmplRatio = repmat(NetPars.AmplRatio, [NetPars.numGroupPerNet,1,1]);

PosiArray = permute(PosiArray, [2,3,1]);
% find the most sensitive neuron
% IdxNeuron = 90;
IdxNeuron = 45;
IdxGroup = [1,3]; % congruent and opposite neurons in network module 1

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(NetPars);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Produce random seeds
seedNoisArray = initRandSeed(NetPars, dimPar, size(parGrid));

%% Bimodal responses
NetStatStruct = struct('OAvgXTime', [], ...
    'OStdXTime', []);
NetStatStruct = repmat(NetStatStruct, size(parGrid));

% fireRate = cell(size(parGrid));
% stdFireRate = cell(size(parGrid));

parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    netpars = parGrid(iterPar);
    netpars.seedNois = seedNoisArray(iterPar);
    
    % Network input
    InputSet = makeNetInput([], netpars);
    
    % Run simulation
    outArgs = struct('InputSet', [], ...
        'NetStat', NetStatStruct(iterPar));
    [InputSet, NetStatStruct(iterPar)] = simDecenNet(InputSet, netpars, outArgs);
end

%% Transform high-dim struct into a scalar struct with high-dim fields
for varName = fieldnames(NetStatStruct)';
    szFields = size(NetStatStruct(1).(varName{1}));
    NetStat = [NetStatStruct.(varName{1})];
    NetStat = reshape(NetStat, [szFields, size(parGrid)]);
    NetStat = permute(NetStat, [ndims(szFields)+1:ndims(NetStat), 1:ndims(szFields)]); % last dim is index of group
    NetEstimRes.(varName{1}) = NetStat;
end
clear varName NetStat NetStatStruct szFields

%% Find the neural weight
% find the neural weight by least square error
Wneu1 = zeros(length(IdxGroup), size(NetPars.AmplRatio, 2)); % neural weight for cue 1
Wneu2 = zeros(length(IdxGroup), size(NetPars.AmplRatio, 2)); % neural weight for cue 2
for iterGroup = 1: length(IdxGroup)
    for iterAmpl = 1: size(NetPars.AmplRatio, 2)
        Y = NetEstimRes.OAvgXTime(iterAmpl, :, 3, IdxNeuron, IdxGroup(iterGroup));
        Y = Y(:);
        X = squeeze(NetEstimRes.OAvgXTime(iterAmpl, :, 1:2, IdxNeuron, IdxGroup(iterGroup)));
        
        X = [X, ones(size(X,1),1)];
        A = (X'*X)\(X'*Y);
        Wneu1(iterGroup, iterAmpl) = A(1);
        Wneu2(iterGroup, iterAmpl) = A(2);
    end
end
clear X Y A

% Rate3E = bsxfun(@times, reshape(Wneu1, 3, 1), fireRate(:,:,1)) + ...
%     bsxfun(@times, reshape(Wneu2, 3, 1), fireRate(:,:,2));
clear iterGroup iterAmpl

%% Save
% if bSave
%     cd D:\CANN_Coding\Data\doubleCANN\singleNeuron
%     str = date;
%     save(['BimodalRes(', str(1:end-5), ').mat'])
% end

%%
% IdxNeuron = 90;
FontSize = 9;

figure;
% cMap = flip(hot(64), 1);
% colormap(cMap)
colormap(jet)

pos = get(gcf, 'position');
pos(4) = pos(3);
pos(2) = pos(2) - 140;
set(gcf, 'position', pos);

cMin = min(reshape(NetEstimRes.OAvgXTime(:, :, :, IdxNeuron, [1,3]), 1,[]));
cMax = max(reshape(NetEstimRes.OAvgXTime(:, :, :, IdxNeuron, [1,3]), 1,[]));

for iterGroup = 1: length(IdxGroup)
    fireRate = squeeze(NetEstimRes.OAvgXTime(:, :, :, IdxNeuron, IdxGroup(iterGroup)));
    
    for iterAmpl = 1: size(NetPars.AmplRatio, 2)
        %         hAxe(i) = subplot(1, nAmpl, i);
        
        hAxe(iterAmpl) = axes('position', [0.15+0.3*(iterAmpl-1), 1.2-0.5*iterGroup, 0.2, 0.2]);
%         hAxe(iterAmpl) = axes('position', [0.1+0.3*(iterAmpl-1), 0.4, 0.2, 0.2]);
        
        fireRatePlot = reshape(squeeze(fireRate(iterAmpl,:,3)), 9,9);
        contourf(PosiArray(:,:,1)', PosiArray(:,:,2)', fireRatePlot); %, ...
        %'levelstep', 4);
        axis xy
%         caxis([min(fireRate(:)), max(fireRate(:))])
        caxis([cMin, cMax])
        set(hAxe(iterAmpl), 'xtick', [], 'ytick',[], 'fontsize', FontSize)
        %     axis square
        daspect([1,1,1])
        %     set(gca,'units','centimeters')
        
        % plot the margical curve on y-axis
        figpos = get(hAxe(iterAmpl), 'position');
        figpos(1) = figpos(1) - 0.07;
        figpos(3) = figpos(3)/4;
        hAxeMar(2*iterAmpl-1) = axes('position', figpos);
        
        Idx = (NetPars.Posi(2,:)==0);
        plot(squeeze(fireRate(iterAmpl, Idx, 1)), PosiArray(1,:,1), ...
            '-ok', 'linew',2, 'markersize', 4, 'markerfacecolor', 'k', 'markeredgecolor', 'none')
        set(gca, 'xlim', [0, max(fireRate(:))], 'ylim', [-180, 180])
        set(gca, 'ytick', -180:90:180,  'xtick', [0, 30]);
        set(gca, 'xaxislocation', 'top', 'color', 'none','fontsize', FontSize, 'linew',1);
        box off
        %     set(gca,'units','centimeters')
        
        % plot the margical curve on x-axis
        figpos = get(hAxe(iterAmpl), 'position');
        figpos(2) = figpos(2) - 0.07;
        figpos(4) = figpos(4)/4;
        hAxeMar(2*iterAmpl) = axes('position', figpos);
        
        Idx = (NetPars.Posi(1,:)==0);
        plot(PosiArray(:,1,2), squeeze(fireRate(iterAmpl, Idx, 2)),  ...
            '-ok', 'linew',2, 'markersize', 4, 'markerfacecolor', 'k', 'markeredgecolor', 'none')
        set(gca, 'ylim', [0, max(fireRate(:))], 'xlim', [-180, 180])
        set(gca, 'xtick', -180:90:180, 'ytick', [0, 30]);
        set(gca, 'color', 'none', 'fontsize', FontSize, 'linew',1)
        box off
        %     set(gca,'units','centimeters')
        
    end
end

% AmplRatio = NetPars.AmplRatio(1,:)./NetPars.AmplRatio(2,:);
% title(['\alpha_1 =' num2str(NetPars.AmplRatio) '\alpha_2'], 'fontsize', FontSize)

for iter = 1: size(NetPars.AmplRatio, 2)
axes(hAxe(iter))    
title(['\alpha_1=' num2str(NetPars.AmplRatio(1,iter)), ...
    ', \alpha_2=' num2str(NetPars.AmplRatio(2,iter))], 'fontsize', FontSize)
end

axes(hAxeMar(1))
ylabel('Cue 1')
axes(hAxeMar(4))
xlabel('Cue 2')

axes(hAxe(3))
hAxe(4) = colorbar;
set(hAxe(3), 'position', [0.15+0.3*(iterAmpl-1), 1.2-0.5*iterGroup, 0.2, 0.2]);
% axes(hAxe(4))
% ylabel('Firing rate', 'fontsize', FontSize)
% set(hAxe(4), 'ytick', 5:10:25, 'fontsize', FontSize)

%%
figure
for i = 1: 2
    hAxe1(i) = subplot(1, 2, i);
    hold on
end
plot(hAxe1(1), NetPars.AmplRatio(1,:)./NetPars.AmplRatio(2,:), Wneu1, 'linew', 2)
plot(hAxe1(1), NetPars.AmplRatio(1,:)./NetPars.AmplRatio(2,:), Wneu2, 'r', 'linew', 2)
plot(hAxe1(2), NetPars.AmplRatio(1,:)./NetPars.AmplRatio(2,:), Wneu1./Wneu2, 'linew', 2)

for i = 1: 2
    axes(hAxe1(i));
    axis square;
end

axes(hAxe1(1))
ylabel('Neural weight', 'fontsize', 9)
% set(gca, 'xlim', [0.1, 0.7], ...
%     'ylim', [0, 1], 'ytick', [0,0.5,1])
xlabel('Cue 1 intensity \alpha_1 (/\alpha_2)', 'fontsize', 9)
legend('w_1','w_2','location', 'southwest')

axes(hAxe1(2))
ylabel('W_1/W_2', 'fontsize', 9)

