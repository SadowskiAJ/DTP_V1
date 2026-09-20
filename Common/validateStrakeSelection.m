function mapped = validateStrakeSelection(strakes,bottoms,tops,bounds,minZ,maxZ,tolerance,retained,detected)
% Place validated internal joints at their nominal boundary indices.
if nargin < 8
    retained = true(size(bounds,1),1);
end
retained = retained(:);
if numel(retained) ~= size(bounds,1) || ...
        any(~isfinite(bounds(retained,:)),'all') || ...
        any(~isnan(bounds(~retained,:)),'all')
    error('validateStrakeSelection:InvalidRetainedBounds', ...
        'Retained joint bounds must be finite; excluded joint bounds must be NaN.')
end
validateattributes(strakes,{'numeric'},{'vector','nonempty','integer','positive','<=',numel(bottoms)})
validateattributes(tolerance,{'numeric'},{'real','scalar','finite','positive'})
strakes = strakes(:);
bottoms = bottoms(:)/1000; tops = tops(:)/1000;
if numel(bottoms) ~= numel(tops) || any(~isfinite([bottoms;tops])) || any(tops <= bottoms) || ...
        size(bounds,2) ~= 2 || any(isinf(bounds(:))) || ...
        any(bounds(:,2) < bounds(:,1)) || numel(unique(strakes)) ~= numel(strakes)
    error('validateStrakeSelection:InvalidBounds','Nominal geometry, selected strakes or fitted joint bounds are invalid.')
end
outside = bottoms(strakes) < minZ | tops(strakes) > maxZ;
if any(outside)
    error('validateStrakeSelection:OutsideRange', ...
        'Selected strakes %s are not wholly inside %.3f--%.3f m. Review cloudStrakes and the height range.', ...
        mat2str(strakes(outside)'),minZ,maxZ)
end
[~,boundaryRows] = selectCloudStrakes(bottoms*1000,tops*1000,minZ,maxZ);
if size(bounds,1) ~= numel(boundaryRows)
    error('validateStrakeSelection:MissingBoundary', ...
        'Expected %g internal joint bounds but received %g. No missing boundary has been inferred.', ...
        numel(boundaryRows),size(bounds,1))
end
% Row S bounds the bottom of nominal segment S. End rows stay missing.
mapped = nan(numel(bottoms)+1,2);
mapped(boundaryRows,:) = bounds;
centres = mean(mapped,2);
if nargin >= 9
    if numel(detected) ~= numel(boundaryRows)
        error('validateStrakeSelection:InvalidDetection','Detected joints must match the nominal boundary rows.')
    end
    % The nominal tolerance applies to detected peaks, not asymmetric ROI midpoints.
    centres(boundaryRows) = detected(:);
end
allowed = min(tolerance,(tops(strakes)-bottoms(strakes))/2);
bad = ~isfinite(centres(strakes)) | ~isfinite(centres(strakes+1)) | ...
    abs(centres(strakes)-bottoms(strakes)) >= allowed | ...
    abs(centres(strakes+1)-tops(strakes)) >= allowed | ...
    mapped(strakes,2) >= mapped(strakes+1,1);
if any(bad)
    error('validateStrakeSelection:AssignmentMismatch', ...
        ['Detected bounds do not correspond to nominal strakes %s. Review Stage 2 and cloudStrakes. ', ...
         'No missing boundary or tower end has been inferred.'],mat2str(strakes(bad)'))
end
end
