function [minZ,maxZ] = resolveHeightRange(inputDir,minZ,maxZ)
% Resolve open height limits from all binary input points, in metres.
% Finite limits are used as supplied. The scan is cached by input signature.
if ~(isnumeric(minZ) && isscalar(minZ) && isreal(minZ) && ...
        isnumeric(maxZ) && isscalar(maxZ) && isreal(maxZ) && ...
        ~isnan(minZ) && ~isnan(maxZ) && minZ < maxZ)
    error('resolveHeightRange:InvalidRange','Height limits must be real scalars with minZ < maxZ.')
end
if isfinite(minZ) && isfinite(maxZ)
    return
end

signature = pointCloudFileSignature(inputDir);
cacheFile = fullfile(inputDir,'S0_PointCloudHeight.mat');
bounds = [];
if isfile(cacheFile)
    cached = load(cacheFile,'inputSignature','heightBounds');
    if isfield(cached,'inputSignature') && isequaln(cached.inputSignature,signature) && ...
            isfield(cached,'heightBounds') && isnumeric(cached.heightBounds) && ...
            isequal(size(cached.heightBounds),[1 2]) && ...
            all(isfinite(cached.heightBounds)) && cached.heightBounds(1) < cached.heightBounds(2)
        bounds = cached.heightBounds;
    end
end
if isempty(bounds)
    bounds = [Inf -Inf];
    for scan = 1:numel(signature)
        fileName = fullfile(inputDir,signature(scan).name);
        fprintf('Reading height range: %s (%g of %g)\n',signature(scan).name,scan,numel(signature))
        scanBounds = binaryHeightBounds(fileName,signature(scan).bytes);
        bounds = [min(bounds(1),scanBounds(1)),max(bounds(2),scanBounds(2))];
    end
    if any(~isfinite(bounds)) || bounds(1) >= bounds(2)
        error('resolveHeightRange:NoHeight','The input points have no finite positive height range.')
    end
    inputSignature = signature;
    heightBounds = bounds;
    save(cacheFile,'inputSignature','heightBounds')
end
if isinf(minZ); minZ = bounds(1); end
if isinf(maxZ); maxZ = bounds(2); end
if minZ >= maxZ
    error('resolveHeightRange:EmptyRange','The resolved height limits do not define a positive range.')
end
fprintf('Height range: %.6f to %.6f m\n',minZ,maxZ)
end

function bounds = binaryHeightBounds(fileName,fileSize)
fileID = fopen(fileName,'rb','ieee-le');
if fileID < 0
    error('resolveHeightRange:FileOpenFailed','Unable to open %s.',fileName)
end
cleanup = onCleanup(@() fclose(fileID));
count = fread(fileID,1,'uint64=>double');
if isempty(count) || count > flintmax || fileSize ~= 8+24*count
    error('resolveHeightRange:InvalidBinary','Invalid point count or file size in %s.',fileName)
end
bounds = [Inf -Inf];
remaining = count;
while remaining > 0
    chunkSize = min(remaining,1e6); % At most 24 MB of coordinates at a time.
    points = fread(fileID,[3 chunkSize],'double=>double');
    if numel(points) ~= 3*chunkSize
        error('resolveHeightRange:ReadFailed','Incomplete coordinate data in %s.',fileName)
    end
    if any(~isfinite(points(:)))
        error('resolveHeightRange:InvalidCoordinates','Nonfinite coordinates in %s.',fileName)
    end
    bounds = [min(bounds(1),min(points(3,:))),max(bounds(2),max(points(3,:)))];
    remaining = remaining-chunkSize;
end
end
