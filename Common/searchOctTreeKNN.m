% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [dists, nearestInds] = searchOctTreeKNN(point, knear, points, ~, ~, ~)
%SEARCHOCTTREEKNN Compatibility wrapper for exact compiled k-NN search.
%
% The original public API is retained, but MATLAB's compiled KD-tree search
% replaces the interpreted octree traversal. Tied neighbours may be returned
% in a different order; distances and neighbour sets are otherwise equivalent.

point = reshape(point,1,[]);
if numel(point) ~= 3 || size(points,2) ~= 3 || isempty(points)
    error('searchOctTreeKNN:InvalidInput', ...
        'point must contain three coordinates and points must be nonempty N-by-3.')
end
if ~(isscalar(knear) && knear >= 1 && knear == floor(knear))
    error('searchOctTreeKNN:InvalidK','knear must be a positive integer.')
end
knear = min(knear,size(points,1));
[nearestInds,dists] = knnsearch(points,point,'K',knear,'Distance','euclidean');
nearestInds = nearestInds(:);
dists = dists(:);
end
