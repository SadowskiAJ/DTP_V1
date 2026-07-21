% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [dists, nearestInds] = searchOctTreeRadius(point, radius, points, ~, ~, ~)
%SEARCHOCTTREERADIUS Compatibility wrapper for exact compiled radius search.

point = reshape(point,1,[]);
if numel(point) ~= 3 || size(points,2) ~= 3 || isempty(points)
    error('searchOctTreeRadius:InvalidInput', ...
        'point must contain three coordinates and points must be nonempty N-by-3.')
end
if ~(isscalar(radius) && isfinite(radius) && radius >= 0)
    error('searchOctTreeRadius:InvalidRadius','radius must be a finite nonnegative scalar.')
end
[indexCells,distanceCells] = rangesearch(points,point,radius,'Distance','euclidean');
nearestInds = indexCells{1}(:);
dists = distanceCells{1}(:);
[dists,order] = sort(dists);
nearestInds = nearestInds(order);
end
