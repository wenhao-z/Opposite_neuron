function seedNoisArray = initRandSeed(NetPars, dimPar, sizeGrid)
% Produce random seeds for scanNetPars.m
% Wen-Hao Zhang, Apr-7, 2016

% sizeGrid: size of the parameter Grid for scanning parameters
for iter = 1: length(dimPar)
    if strcmp(dimPar(iter).namePar, 'cueCond')
        dimCueCond = iter;
        break
    end
end

switch NetPars.flagSeed
    case 'sameSeed'
        seedNoisArray = repmat(NetPars.seedNois, sizeGrid);
    case 'SameSeedCueCond'
        sizeGrid(dimCueCond) = 1;
        sizeGridRep = ones(size(sizeGrid));
        sizeGridRep(dimCueCond) = 3;
        
        seedNoisArray = repmat(rand(sizeGrid)*1e5, sizeGridRep);
    case 'diffSeed'
        sizeGrid(dimCueCond) = 3;
        seedNoisArray = rand(sizeGrid)*1e5;
end



% Followings are old version which generate the seeds for external and
% internal noises respectively.

% switch NetPars.flagSeed
%     case 'sameSeed'
%         seedNoisArray = repmat(NetPars.seedIntNois, sizeGrid);
%         seedExtNoisArray = repmat(NetPars.seedExtNois, sizeGrid);
%     case 'SameSeedCueCond'
%         sizeGrid(dimCueCond) = 1;
%         sizeGridRep = ones(size(sizeGrid));
%         sizeGridRep(dimCueCond) = 3;
%         
%         seedNoisArray = repmat(rand(sizeGrid)*1e5, sizeGridRep);
%         seedExtNoisArray = repmat(rand(sizeGrid)*1e5, sizeGridRep);
%     case 'diffSeed'
%         sizeGrid(dimCueCond) = 3;
%         seedNoisArray = rand(sizeGrid)*1e5;
%         seedExtNoisArray = rand(sizeGrid)*1e5;
% end

