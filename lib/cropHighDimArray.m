function TestRes = cropHighDimArray(TestRes, namePar, IdxPars, dimOrder)
% Crop the high-dim array on each dim of parameters
% INPUT:
% TestRes: high-dim arry to be cropped
% namePar: the name of each dimension of TestRes
% IdxPars: the index of cropping on each dimension of TestRes
% dimOrder: the first several dimensions of TestRes after cropping

% Wen-Hao Zhang, June-18,2017
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

for parName = fieldnames(IdxPars)';
    if strcmp(parName{1}, 'IdxNeuronGroup')
        % Caution: make sure the LAST dim is index of group (check with the code scanNetPars.)
        IdxDim = length(namePar) + 1;
    else
        IdxDim = cellfun(@(x) strcmp(x, parName{1}), namePar);
        IdxDim = find(IdxDim);
    end
    if isempty(IdxDim)
        continue;
    end
    
    TestRes = permute(TestRes, [IdxDim, setdiff(1:ndims(TestRes), IdxDim)]);
    szTestRes_Perm = size(TestRes);
    TestRes = reshape(TestRes, size(TestRes,1), []);
    TestRes = TestRes(IdxPars.(parName{1}),:); % Crop on this dim
    TestRes = reshape(TestRes, [size(TestRes,1), szTestRes_Perm(2:end)]);
    TestRes = ipermute(TestRes, [IdxDim, setdiff(1:ndims(TestRes), IdxDim)]);
end

% Permute the dimension of the cropped high-dim array as indicated by dimOrder
TestRes = permute(TestRes, [dimOrder, setdiff(1:ndims(TestRes), dimOrder)]);
switch length(dimOrder)
    case 1
        TestRes = reshape(TestRes, size(TestRes, 1), []);
    otherwise
        TestRes = reshape(TestRes, size(TestRes, 1), size(TestRes, 2), []);
end
end