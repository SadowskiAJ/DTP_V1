% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function compatible = isCompatibleCache(fileName,expectedConfig)
%ISCOMPATIBLECACHE True when a MAT-file contains the expected cache metadata.

compatible = false;
if ~isfile(fileName)
    return
end
variables = whos('-file',fileName);
if ~any(strcmp({variables.name},'cacheConfig'))
    return
end
stored = load(fileName,'cacheConfig');
compatible = isequaln(stored.cacheConfig,expectedConfig);
end
