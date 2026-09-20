% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [THETg, Zg, RBIG, failure] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd, plotting)
%S4_GRIDDEDOUTERSHELLSURFACERECONSTRUCTION Reconstruct a regular cylindrical surface.
% An unusable cone returns empty grids and a failure reason for the caller.
THETg = []; Zg = []; RBIG = []; failure = [];

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
% Points are already aligned to the tower axis. Their centroid is biased
% towards the observed side when circumferential coverage is uneven.
initial = [0,0,initialRadius,initialRadius];
options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5e5, 'FunctionTolerance', 1e-14, ...
    'StepTolerance', 1e-14, 'MaxIterations', 5000);
[coneParameters,~,~,exitFlag] = lsqnonlin(@(b) basicCone(b,points), initial, [], [], options);
if exitFlag <= 0 || any(~isfinite(coneParameters))
    failure = struct('identifier','S4_GriddedOuterShellSurfaceReconstruction:ConeFitFailed', ...
        'reason',sprintf('Cone fit did not converge (exit flag %g).',exitFlag));
    if nargout < 4; warning(failure.identifier,'Strake %g: %s',S,failure.reason); end
    return
end

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
if any(~isfinite(rloc) | rloc <= 0)
    failure = struct('identifier','S4_GriddedOuterShellSurfaceReconstruction:InvalidCone', ...
        'reason','The fitted cone produced a nonfinite or nonpositive radius.');
    if nargout < 4; warning(failure.identifier,'Strake %g: %s',S,failure.reason); end
    return
end

nTheta = max(3,ceil(2*pi*min(rloc)/usedSpacing));
theta = (0:nTheta-1)*(2*pi/nTheta);
[THETg,Zg] = meshgrid(theta,z);
RBIG = rloc(:).*ones(1,nTheta);
XNEW = RBIG.*cos(THETg)+coneParameters(1);
YNEW = RBIG.*sin(THETg)+coneParameters(2);

searchRadius = usedSpacing*searchRadiusGridSpacingSF;
queries = [XNEW(:),YNEW(:),Zg(:)];
gridRadii = inverseDistanceValues(searcher,queries,RR,searchRadius,p);
RBIG = reshape(gridRadii,size(XNEW));
fprintf('Stage 4 strake %g: %.1f%% of grid nodes unsupported (search radius %.2f mm).\n', ...
    S,100*nnz(isnan(RBIG))/numel(RBIG),searchRadius)

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

function values = inverseDistanceValues(searcher,queries,sourceValues,radius,p)
values = nan(size(queries,1),1);
batchSize = 5000;
for first = 1:batchSize:size(queries,1)
    batch = first:min(first+batchSize-1,size(queries,1));
    % A fixed local support radius leaves unmeasured areas as NaN.
    [indices,distances] = rangesearch(searcher,queries(batch,:),radius);
    for cellIndex = find(~cellfun('isempty',indices))'
        sourceIndices = indices{cellIndex};
        sourceDistances = distances{cellIndex};
        zeroDistance = sourceDistances <= 64*eps(max(1,max(sourceDistances)));
        if any(zeroDistance)
            values(batch(cellIndex)) = mean(sourceValues(sourceIndices(zeroDistance)));
        else
            weights = 1./sourceDistances(:).^p;
            values(batch(cellIndex)) = ...
                sum(sourceValues(sourceIndices(:)).*weights)/sum(weights);
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
% Normalise by observed support so NaNs do not contaminate nearby data.
% Restore the original mask afterwards: smoothing must not fill gaps.
support = isfinite(extended);
extended(~support) = 0;
extended = imgaussfilt(extended,sigma,'FilterSize',filterSize,'Padding','replicate', ...
    'FilterDomain','spatial');
weights = imgaussfilt(double(support),sigma,'FilterSize',filterSize,'Padding','replicate', ...
    'FilterDomain','spatial');
extended = extended./weights;
filtered = extended(:,halfWidth+(1:nColumns));
filtered(~isfinite(values)) = NaN;
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
