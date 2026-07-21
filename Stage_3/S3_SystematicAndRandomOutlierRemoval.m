% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [TT_REM, ZZ_REM, RR_REM, TT_OUT_BELOW, ZZ_OUT_BELOW, RR_OUT_BELOW, TT_OUT_ABOVE, ZZ_OUT_ABOVE, RR_OUT_ABOVE, DISTPQ] = S3_SystematicAndRandomOutlierRemoval(TT, ZZ, RR, RMEAN, tbinSize, zbinSize, cutOff, plotting)
%S3_SYSTEMATICANDRANDOMOUTLIERREMOVAL Remove points far from local PCA planes.

if nargin < 8
    plotting = false;
end
TT = TT(:);
ZZ = ZZ(:);
RR = RR(:);
if numel(TT) ~= numel(ZZ) || numel(TT) ~= numel(RR) || isempty(TT)
    error('S3_SystematicAndRandomOutlierRemoval:InvalidInput', ...
        'TT, ZZ and RR must be nonempty vectors of equal length.')
end
if any(~isfinite([TT; ZZ; RR])) || RMEAN <= 0 || ...
        tbinSize <= 0 || zbinSize <= 0 || cutOff < 0
    error('S3_SystematicAndRandomOutlierRemoval:InvalidParameter', ...
        'Coordinates and window parameters must be finite and physically valid.')
end

originalTheta = mod(TT, 2*pi);
originalZ = ZZ;
originalR = RR;
nOriginal = numel(TT);
minZZ = min(ZZ);
ZZ = ZZ-minZZ;

% Scale theta to arc length and extend both periodic boundaries using only
% the original points. Computing both selections before either append avoids
% duplicating the first extension back into the physical domain.
circumference = 2*pi*RMEAN;
TT = originalTheta*RMEAN;
rightOriginalInds = find(TT > circumference-tbinSize);
leftOriginalInds = find(TT < tbinSize);

rightAddedInds = nOriginal+(1:numel(rightOriginalInds));
leftAddedInds = nOriginal+numel(rightOriginalInds)+(1:numel(leftOriginalInds));
TT = [TT; TT(rightOriginalInds)-circumference; TT(leftOriginalInds)+circumference];
ZZ = [ZZ; ZZ(rightOriginalInds); ZZ(leftOriginalInds)];
RR = [RR; RR(rightOriginalInds); RR(leftOriginalInds)];

thetaStarts = windowStarts(min(TT), max(TT), tbinSize);
zStarts = windowStarts(min(ZZ), max(ZZ), zbinSize);
[thetaStartGrid,zStartGrid] = meshgrid(thetaStarts,zStarts);
windowCentres = [thetaStartGrid(:)+tbinSize/2, zStartGrid(:)+zbinSize/2];

% A compiled KD-tree supplies a conservative circular candidate set for
% each rectangular window; the exact rectangular condition is then applied.
scaledPoints = [TT/tbinSize, ZZ/zbinSize];
scaledCentres = [windowCentres(:,1)/tbinSize, windowCentres(:,2)/zbinSize];
searcher = createns(scaledPoints, 'NSMethod', 'kdtree', 'Distance', 'euclidean');
candidateIndices = rangesearch(searcher, scaledCentres, sqrt(0.5));

DISTPQ = zeros(size(TT));
inadequateWindows = 0;
for windowIndex = 1:size(windowCentres,1)
    inds = candidateIndices{windowIndex};
    if isempty(inds)
        inadequateWindows = inadequateWindows+1;
        continue
    end
    centre = windowCentres(windowIndex,:);
    inside = abs(TT(inds)-centre(1)) < tbinSize/2 & ...
        abs(ZZ(inds)-centre(2)) < zbinSize/2;
    inds = inds(inside);
    if numel(inds) <= 6
        inadequateWindows = inadequateWindows+1;
        continue
    end

    coordinates = [TT(inds), ZZ(inds), RR(inds)];
    coefficients = pca(coordinates);
    normal = coefficients(:,3);
    if normal(3) < 0
        normal = -normal;
    end

    initialDistances = (coordinates-coordinates(1,:))*normal;
    [~,medianPointIndex] = min(abs(initialDistances-median(initialDistances)));
    distances = (coordinates-coordinates(medianPointIndex,:))*normal;

    existing = DISTPQ(inds);
    replace = abs(distances) > abs(existing);
    existing(replace) = distances(replace);
    DISTPQ(inds) = existing;
end
if inadequateWindows > 0
    warning('S3_SystematicAndRandomOutlierRemoval:SparseWindows', ...
        '%g sliding windows contained six or fewer points.', inadequateWindows)
end

% Transfer the maximum distance observed in each periodic copy back to its
% original point, then discard all extension points by direct indexing.
originalDistances = DISTPQ(1:nOriginal);
originalDistances = mergePeriodicDistances(originalDistances, DISTPQ, ...
    rightOriginalInds, rightAddedInds);
originalDistances = mergePeriodicDistances(originalDistances, DISTPQ, ...
    leftOriginalInds, leftAddedInds);
DISTPQ = originalDistances;
TT = originalTheta;
ZZ = originalZ;
RR = originalR;

inlier = abs(DISTPQ) <= cutOff;
above = DISTPQ > cutOff;
below = DISTPQ < -cutOff;
TT_REM = TT(inlier);
ZZ_REM = ZZ(inlier);
RR_REM = RR(inlier);
TT_OUT_ABOVE = TT(above);
ZZ_OUT_ABOVE = ZZ(above);
RR_OUT_ABOVE = RR(above);
TT_OUT_BELOW = TT(below);
ZZ_OUT_BELOW = ZZ(below);
RR_OUT_BELOW = RR(below);

if plotting
    figure
    hold on
    scatter3(TT_REM, ZZ_REM, RR_REM, 20, '.k')
    scatter3(TT_OUT_BELOW, ZZ_OUT_BELOW, RR_OUT_BELOW, 20, '.r')
    scatter3(TT_OUT_ABOVE, ZZ_OUT_ABOVE, RR_OUT_ABOVE, 20, '.b')
    grid on
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

function starts = windowStarts(minimumValue, maximumValue, width)
available = maximumValue-minimumValue;
if available <= width
    starts = minimumValue;
    return
end
starts = minimumValue:width/3:(maximumValue-width);
if starts(end) < maximumValue-width-64*eps(maximumValue)
    starts(end+1) = maximumValue-width;
end
end

function original = mergePeriodicDistances(original, allDistances, originalInds, addedInds)
if isempty(originalInds)
    return
end
replace = abs(allDistances(addedInds)) > abs(original(originalInds));
original(originalInds(replace)) = allDistances(addedInds(replace));
end
