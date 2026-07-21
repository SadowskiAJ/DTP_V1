% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function convertPTS2BIN(tower)
%CONVERTPTS2BIN Convert every missing or invalid PTS cloud to binary form.

commonDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(commonDir);
inputDir = fullfile(repoRoot, ['Input_', char(tower)]);
converter = fullfile(commonDir, 'PTS2BIN.exe');

if ~isfolder(inputDir)
    error('convertPTS2BIN:InputDirectoryMissing', ...
        'Input directory does not exist: %s', inputDir)
end
if ~isfile(converter)
    error('convertPTS2BIN:ConverterMissing', ...
        'Converter executable does not exist: %s', converter)
end

ptsFiles = dir(fullfile(inputDir, '*.pts'));
if isempty(ptsFiles)
    if isempty(dir(fullfile(inputDir, '*.bin')))
        error('convertPTS2BIN:NoInputFiles', ...
            'No .pts or .bin files were found in %s.', inputDir)
    end
    disp('No .pts files found; using existing .bin files.')
    return
end

for j = 1:numel(ptsFiles)
    ptsPath = fullfile(ptsFiles(j).folder, ptsFiles(j).name);
    [~, stem] = fileparts(ptsPath);
    fileRoot = fullfile(ptsFiles(j).folder, stem);
    binPath = [fileRoot, '.bin'];

    if isValidBin(binPath)
        fprintf('Valid .bin file %g of %g already exists: %s\n', ...
            j, numel(ptsFiles), [stem, '.bin'])
        continue
    end

    fprintf('Converting %g of %g: %s\n', j, numel(ptsFiles), ptsFiles(j).name)
    command = sprintf('"%s" "%s"', converter, fileRoot);
    [status, output] = system(command);
    if status ~= 0
        error('convertPTS2BIN:ConversionFailed', ...
            'PTS2BIN failed for %s (exit status %d):\n%s', ptsPath, status, output)
    end
    if ~isValidBin(binPath)
        error('convertPTS2BIN:InvalidOutput', ...
            'PTS2BIN reported success but produced an invalid file: %s', binPath)
    end
end
end

function valid = isValidBin(fileName)
valid = false;
info = dir(fileName);
if isempty(info) || info.bytes < 8
    return
end

fileID = fopen(fileName, 'rb', 'ieee-le');
if fileID < 0
    return
end
cleaner = onCleanup(@() fclose(fileID));
nPoints = fread(fileID, 1, 'uint64=>double');
valid = isscalar(nPoints) && isfinite(nPoints) && nPoints >= 0 && ...
    info.bytes == 8 + 3*8*nPoints;
end
