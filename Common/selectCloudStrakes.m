function [strakes,boundaryRows] = selectCloudStrakes(bottoms,tops,minZ,maxZ)
% Select complete segments bounded by two internal nominal junctions.
bottoms = bottoms(:)/1000; tops = tops(:)/1000;
if numel(bottoms) ~= numel(tops) || isempty(bottoms) || ...
        any(~isfinite([bottoms;tops])) || any(tops <= bottoms) || ...
        any(diff(bottoms) <= 0) || any(tops(1:end-1) ~= bottoms(2:end))
    error('selectCloudStrakes:InvalidGeometry', ...
        'Nominal segments must be finite, ordered and share consecutive boundaries.')
end
validateattributes(minZ,{'numeric'},{'scalar','real','finite'})
validateattributes(maxZ,{'numeric'},{'scalar','real','finite','>',minZ})
boundaryRows = (2:numel(bottoms))';
boundaryRows = boundaryRows(bottoms(boundaryRows) > minZ & bottoms(boundaryRows) < maxZ);
available = false(numel(bottoms)+1,1);
available(boundaryRows) = true;
strakes = find(available(1:end-1) & available(2:end))';
if isempty(strakes)
    error('selectCloudStrakes:NoCompleteStrakes', ...
        'No strake has two internal junctions inside %.3f--%.3f m. Extend the height range.',minZ,maxZ)
end
end
