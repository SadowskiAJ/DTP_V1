% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function signature = pointCloudFileSignature(inputDir)
%POINTCLOUDFILESIGNATURE Metadata used to detect stale cached stages.

files = dir(fullfile(inputDir,'*.bin'));
if isempty(files)
    error('pointCloudFileSignature:NoBinaryClouds', ...
        'No .bin point clouds were found in %s.',inputDir)
end
[~,order] = sort({files.name});
files = files(order);
signature = struct('name',{files.name},'bytes',num2cell([files.bytes]), ...
    'modified',{files.date});
end
