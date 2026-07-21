% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [X, Y, Z, T, R] = postRegDataLoader(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower)
%POSTREGDATALOADER Load, downsample and apply the Stage 1 transformations.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(repoRoot, ['Input_', tower]);
files = dir(fullfile(inputDir, '*.bin'));
nClouds = numel(files);
if nClouds == 0
    error('postRegDataLoader:NoClouds', 'No .bin clouds were found in %s.', inputDir)
end

clouds = S1_CLOUD(fullfile(files(1).folder, files(1).name), 1);
for j = 2:nClouds
    clouds(j) = S1_CLOUD(fullfile(files(j).folder, files(j).name), j); %#ok<AGROW>
end

Xparts = cell(1,nClouds);
Yparts = cell(1,nClouds);
Zparts = cell(1,nClouds);
for j = 1:nClouds
    clouds(j) = clouds(j).CLOUDDOWNRND(nperc);
    [clouds(j).X, clouds(j).Y, clouds(j).Z, clouds(j).T, clouds(j).R] = ...
        bestFitConeTransform(clouds(j).X, clouds(j).Y, clouds(j).Z, ...
        minZ, bestFitParamsPreReg, false);
    clouds(j) = clouds(j).CLOUDMERGEITER(regParams, iUse);
    Xparts{j} = clouds(j).XNEW(:);
    Yparts{j} = clouds(j).YNEW(:);
    Zparts{j} = clouds(j).ZNEW(:);
end

X = vertcat(Xparts{:});
Y = vertcat(Yparts{:});
Z = vertcat(Zparts{:});
selection = Z >= minZ & Z <= maxZ;
X = X(selection);
Y = Y(selection);
Z = Z(selection);

T = mod(atan2(Y,X), 2*pi);
R = hypot(X,Y);
end
