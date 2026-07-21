% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [THETg, Zg, RBIG] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd, plotting)
%S4_GRIDDEDOUTERSHELLSURFACERECONSTRUCTION Reconstruct a regular cylindrical surface.

if nargin < 15
    plotting = false;
end
TT = TT_REM(:);
RR = RR_REM(:);
ZZ = ZZ_REM(:)-botExtent;
if numel(TT) ~= numel(RR) || numel(TT) ~= numel(ZZ) || numel(TT) < 3
    error('S4_GriddedOuterShellSurfaceReconstruction:InvalidInput', ...
        'At least three theta, radius and elevation coordinate triples are required.')
end
if any(~isfinite([TT; RR; ZZ])) || topExtent <= botExtent || ...
        searchRadiusGridSpacingSF <= 0 || p <= 0 || smoothingFilterStd < 0
    error('S4_GriddedOuterShellSurfaceReconstruction:InvalidParameter', ...
        'Coordinates and reconstruction parameters must be finite and physically valid.')
end

XX = RR.*cos(TT);
YY = RR.*sin(TT);
points = [XX,YY,ZZ];
searcher = createns(points, 'NSMethod', 'kdtree', 'Distance', 'euclidean');

% Characteristic point spacing from exact batched 1-NN queries.
nRandom = min(1000,numel(TT));
randomIndices = randperm(numel(TT),nRandom);
[~,nearestDistances] = knnsearch(searcher, points(randomIndices,:), 'K', 2);
nearestDistances = nearestDistances(:,2);
nearestDistances = nearestDistances(isfinite(nearestDistances) & nearestDistances > 0);
if isempty(nearestDistances)
    error('S4_GriddedOuterShellSurfaceReconstruction:DegenerateSpacing', ...
        'A positive characteristic point spacing could not be determined.')
end
usedSpacing = prctile(nearestDistances,95);

% Best-fitting unrotated cone used to centre the nominal reconstruction grid.
initialRadius = median(RR);
initial = [mean(XX),mean(YY),initialRadius,initialRadius];
options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5e5, 'FunctionTolerance', 1e-14, ...
    'StepTolerance', 1e-14, 'MaxIterations', 5000);
coneParameters = lsqnonlin(@(b) basicCone(b,points), initial, [], [], options);

height = topExtent-botExtent;
nZ = max(2,ceil(height/usedSpacing)+1);
z = linspace(0,height,nZ);
zOne = min(ZZ);
zTwo = max(ZZ);
if zTwo <= zOne
    error('S4_GriddedOuterShellSurfaceReconstruction:ZeroHeight', ...
        'The retained can points have no positive elevation range.')
end
slope = (coneParameters(4)-coneParameters(3))/(zTwo-zOne);
intercept = coneParameters(3)-slope*zOne;
rloc = slope*z+intercept;
if any(rloc <= 0)
    error('S4_GriddedOuterShellSurfaceReconstruction:InvalidCone', ...
        'The fitted cone produced a nonpositive radius.')
end

nTheta = max(3,ceil(2*pi*min(rloc)/usedSpacing));
theta = (0:nTheta-1)*(2*pi/nTheta);
[THETg,Zg] = meshgrid(theta,z);
RBIG = rloc(:).*ones(1,nTheta);
XNEW = RBIG.*cos(THETg)+coneParameters(1);
YNEW = RBIG.*sin(THETg)+coneParameters(2);

searchRadius = usedSpacing*searchRadiusGridSpacingSF;
queries = [XNEW(:),YNEW(:),Zg(:)];
gridRadii = inverseDistanceValues(searcher,queries,RR,searchRadius,p,points);
RBIG = reshape(gridRadii,size(XNEW));

if plotting
    figure
    hold on
    scatter3(mod(TT,2*pi),ZZ,RR,1,'.k')
    surf(THETg,Zg,RBIG)
    grid on
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    title(sprintf('Surface reconstruction of strake at %g--%g m, before smoothing', ...
        Z_BOTS(S)/1000,Z_TOPS(S)/1000),'interpreter','latex')
end

bendingWaveLength = 2.444*sqrt(mean([R0_BOTS(S),R0_TOPS(S)])*THICKS(S));
filterSize = floor(bendingWaveLength/5/usedSpacing);
filterSize = max(3,2*floor(filterSize/2)+1);
RBIG = periodicGaussianFilter(RBIG,smoothingFilterStd,filterSize);

if plotting
    figure
    surf(THETg,Zg,RBIG)
    grid on
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    title(sprintf('Surface reconstruction of strake at %g--%g m, after smoothing', ...
        Z_BOTS(S)/1000,Z_TOPS(S)/1000),'interpreter','latex')
end
end

function values = inverseDistanceValues(searcher,queries,sourceValues,initialRadius,p,points)
values = nan(size(queries,1),1);
batchSize = 5000;
maximumRadius = 2*norm(max(points,[],1)-min(points,[],1));
for first = 1:batchSize:size(queries,1)
    batch = first:min(first+batchSize-1,size(queries,1));
    unresolved = (1:numel(batch))';
    radius = initialRadius;
    while ~isempty(unresolved)
        [indices,distances] = rangesearch(searcher,queries(batch(unresolved),:),radius);
        found = ~cellfun('isempty',indices);
        foundCellIndices = find(found);
        for localIndex = 1:numel(foundCellIndices)
            cellIndex = foundCellIndices(localIndex);
            resultPosition = unresolved(cellIndex);
            sourceIndices = indices{cellIndex};
            sourceDistances = distances{cellIndex};
            zeroDistance = sourceDistances <= 64*eps(max(1,max(sourceDistances)));
            if any(zeroDistance)
                values(batch(resultPosition)) = mean(sourceValues(sourceIndices(zeroDistance)));
            else
                weights = 1./sourceDistances(:).^p;
                values(batch(resultPosition)) = ...
                    sum(sourceValues(sourceIndices(:)).*weights)/sum(weights);
            end
        end
        unresolved = unresolved(~found);
        if isempty(unresolved)
            break
        end
        radius = 2*radius;
        if ~isfinite(radius) || radius > maximumRadius
            error('S4_GriddedOuterShellSurfaceReconstruction:NeighbourSearchFailed', ...
                'Unable to find reconstruction neighbours within the point-cloud domain.')
        end
    end
end
end

function filtered = periodicGaussianFilter(values,sigma,filterSize)
if sigma == 0
    filtered = values;
    return
end
halfWidth = (filterSize-1)/2;
nColumns = size(values,2);
leftIndices = mod((-halfWidth:-1),nColumns)+1;
rightIndices = mod((0:halfWidth-1),nColumns)+1;
extended = [values(:,leftIndices),values,values(:,rightIndices)];
extended = imgaussfilt(extended,sigma,'FilterSize',filterSize,'Padding','replicate');
filtered = extended(:,halfWidth+(1:nColumns));
end

function perpImp = basicCone(b,x)
X = x(:,1);
Y = x(:,2);
Z = x(:,3)-min(x(:,3));
height = max(Z);
targetRadius = (b(3)-b(4))/height*(height-Z)+b(4);
beta = atan2(b(3)-b(4),height);
radialImp = targetRadius-hypot(X-b(1),Y-b(2));
perpImp = radialImp.*cos(beta);
end
