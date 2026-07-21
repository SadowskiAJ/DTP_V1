% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Code measures the circumferential dimple tolerance. 

% Measurements are to be taken with the gauge placed on the outside of the
% shell only. There is an option to use the gauge on the inside of the
% shell, but is not accurate as the gauge will always pass through the
% shell wall in such situations.

% An additional check is provided via the selfIntersectionCheck function
% which makes sure no part of the curved gauge (especially from gaugeStart
% to snapOne) intersects with the shell wall after gauge rotation. If an
% intersection is found then the gauge positioning is not feasible and a
% tolerance is not measured. This additional check has been coded for an
% outside shell gauge placement. For nominally circular cross sections with
% small imperfections, this check is unlikely to be needed. This check does
% NOT identify any infeasible gauge placements for the sample WTST dataset.

clear;
clc;
close all

%% Nominal geometry
tower = 'ReleaseTower'; % Tower to be processed
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(fullfile(repoRoot,'Common'))
inputDir = fullfile(repoRoot,['Input_',tower]);
load(fullfile(inputDir,['S6_Mesh_',tower,'.mat']))

%% Extracting joint locations
STRAKES = find(N2_MESH ~= 0); % Finding strakes which have a meshed surface
STRAKES = STRAKES(:).'; % Use a row vector irrespective of N2_MESH orientation
N2_MESH = N2_MESH(STRAKES); % Extracting number of elements in available strakes
N2Cum = cumsum(N2_MESH); % Performing cumulative sum to extract values from overall meshed surface

%% Extracting nominal shell
% req = @(z, S) ((R0_TOPS(S) - R0_BOTS(S))/(Z_TOPS(S) - Z_BOTS(S))*(z - Z_BOTS(S)) + R0_BOTS(S)); % Equation defining nominal radius of shell with nominal z coordinates of joints
req2 = @(z, S, zbot, ztop) ((R0_TOPS(S) - R0_BOTS(S))/(ztop - zbot)*(z - zbot) + R0_BOTS(S)); % Equation defining nominal radius of shell with identified z coordinates of joints (more accurate)
z = Z_MESH(:,1); % Z coordinates
theta = THET_MESH(1,:); % THET coordinates
rnom = zeros(size(R_MESH,1),1); % Nominal radial coordinates
tnom = rnom; % Nominal thickness
lsnom = rnom; % Nominal shell-segment length used by the long gauge
startInd = 1;
for j = 1:length(STRAKES)
    S = STRAKES(j);
    inds = startInd:N2Cum(j); % Indicies of a particular strake
    % rnom(inds) = req(z(inds), S); % Less accurate
    rnom(inds) = req2(z(inds), S, z(inds(1)), z(inds(end)));
    tnom(inds) = THICKS(S)*ones(size(inds));
    % EN 1993-1-6 defines ls between shell-segment boundaries. The nominal
    % strake length is the released default; replace this assignment with
    % design segment lengths where a shell segment spans multiple strakes.
    lsnom(inds) = (Z_TOPS(S)-Z_BOTS(S))*ones(size(inds));
    startInd = N2Cum(j)+1;
end

rnom = rnom - tnom/2; % Correcting to obtain nominal shell midsurface
radialScale = (R_MESH-tnom/2)./R_MESH;
X_MESH = X_MESH.*radialScale; % Converting reconstructed outer surface to midsurface
Y_MESH = Y_MESH.*radialScale;
R_MESH = R_MESH-tnom/2;

%% Entire surface evaluation
plotting = false; % Whether plots are generated or not. Turn off for looped runs!!!

shortGaugeLength = 4*sqrt(rnom.*tnom);
longGaugeLength = max(rnom,2.3*(lsnom.^2.*rnom.*tnom).^(1/4));
if any(shortGaugeLength > pi*rnom) || any(longGaugeLength > pi*rnom)
    error('circumferentialDimpleTolerance:GaugeExceedsSemicircle', ...
        ['Circumferential gauge length must not exceed a semicircle. ', ...
         'Check the nominal radius, thickness and shell-segment length ls.'])
end

fprintf('Evaluating short-wave circumferential dimple tolerance U0x.\n')
U0x = evaluateCircumferentialTolerance(X_MESH,Y_MESH,rnom, ...
    shortGaugeLength,plotting);
fprintf('Evaluating long-wave circumferential dimple tolerance U0theta.\n')
U0theta = evaluateCircumferentialTolerance(X_MESH,Y_MESH,rnom, ...
    longGaugeLength,plotting);

%% Plotting short-wave tolerance variation of U0x
shortFigure = plotToleranceHeatmap(THET_MESH,Z_MESH,R_MESH,U0x, ...
    'Short-wave circumferential dimple tolerance $U_{0x}$', ...
    '$U_{0x}$ [-]');
saveToleranceAssessmentFigure(shortFigure,inputDir, ...
    ['Tolerance_CircumferentialDimple_U0x_',tower])

%% Plotting long-wave tolerance variation of U0theta
longFigure = plotToleranceHeatmap(THET_MESH,Z_MESH,R_MESH,U0theta, ...
    'Long-wave circumferential dimple tolerance $U_{0\theta}$', ...
    '$U_{0\theta}$ [-]');
saveToleranceAssessmentFigure(longFigure,inputDir, ...
    ['Tolerance_CircumferentialDimple_U0theta_',tower])

%% FUNCTIONS
function U0 = evaluateCircumferentialTolerance(X_MESH,Y_MESH,rnom, ...
    gaugeLength,plotting)
U0 = nan(size(X_MESH)); % Unassessed positions remain explicitly undefined
parfor k = 1:size(X_MESH,1)
    fprintf('Examining k = %g of %g\n',k,size(X_MESH,1))
    xv = X_MESH(k,:); yv = Y_MESH(k,:);
    U0row = nan(1,size(X_MESH,2));

    % Best fitting circle to ensure circular cross section has the most
    % accurate centre
    [xv, yv, thetv] = bestFit(xv, yv);
    xv = [xv xv]; yv = [yv yv]; thetv = [thetv thetv+2*pi];
    for j = 1:size(X_MESH,2)

        % Tolerance measurement for a given start point as described in
        % pseudocode
        [lgx, lgxtheta, x1, y1, xcent, ycent, ~, snapOneInd, searchLimit] = defineCircGauge(xv, yv, thetv, rnom, j, k, gaugeLength(k));
        while true
            [snapTwoPos, interPairs, ~, snapOneX, snapOneY, gaugePointX, gaugePointY, gaugeEndX, gaugeEndY, snapTwoPosX, snapTwoPosY, ~] = computeShellGaugeCorrespondences(xv, yv, x1, y1, xcent, ycent, lgxtheta, rnom, searchLimit, snapOneInd, k, plotting);
            [x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, tolMeasureRange, gaugeSnapInd, snapOneInd, snapOneIndOld, endReached] = applyGaugeRotations(gaugePointX, gaugePointY, snapTwoPosX, snapTwoPosY, snapOneInd, snapOneX, snapOneY, snapTwoPos, x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, xv, yv, interPairs, plotting);
            [measureTol] = selfIntersectionCheck(x1, y1, xcent, ycent, snapOneX, snapOneY, xv, yv, rnom, j, k, snapOneIndOld, lgxtheta, plotting);
            if measureTol
                [U0Tol] = measureTolerance(x1, y1, xcent, ycent, snapOneX, snapOneY, gaugeEndX, gaugeEndY, gaugeSnapInd, tolMeasureRange, xv, yv, interPairs, rnom, k, lgx, lgxtheta, snapTwoPos, snapOneIndOld, plotting);
                tolMeasureRange(tolMeasureRange>size(X_MESH,2)) = tolMeasureRange(tolMeasureRange>size(X_MESH,2)) - size(X_MESH,2); % Allowing tolerance measurements to wrap around to start of cross section
                candidate = abs(U0Tol(:));
                current = U0row(tolMeasureRange).';
                replace = isfinite(candidate) & ...
                    (~isfinite(current) | candidate > current);
                current(replace) = candidate(replace);
                U0row(tolMeasureRange) = current.';
            end

            if endReached
                break
            end
        end

    end
    U0(k,:) = U0row;
end
end

function [lgx, lgxtheta, x1, y1, xcent, ycent, perVec, snapOneInd, searchLimit] = defineCircGauge(xv, yv, thetv, rnom, j, k, lgx)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Extract coordinates of gaugeStart from shellPoints
% Initialise snapOneInd, to j representing the start of the gauge index. 
% Compute the gauge length lgx or lgθ using the nominal properties.
% Compute angle swept by the gauge θl.
% Compute normalised vector from gaugeStart to nominal cross-section centre (0,0)
% Compute coordinates of gaugeCentre assuming it lies on vector from gaugeStart to (0,0)
% Compute searchLimit, the index of the point π radians away from snapOneInd
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

r = rnom(k);

x1 = xv(j); y1 = yv(j); % Start of gauge
snapOneInd = j;

lgxtheta = lgx/r; % Angle swept by rod

% Centroid of circle
perVec = [x1; y1]; % Perpendicular from (x1, y1)
perVec = perVec/norm(perVec); % Normalised perpendicular from (x1, y1)

xcent = x1 - r*perVec(1); % X Location of centre of arc
ycent = y1 - r*perVec(2); % X Location of centre of arc

% Not all points along the circumference are considered for snapping. Using
% points further than pi away from snapOne in an anticlockwise direction
% can cause the gauge to potentially rotate the wrong way and snap onto the
% point immediately preceeding the start point (depending on the magnitude
% of imperfections) which is incorrect. Thus points more than pi away from
% snapOne are not to be considered in the algorithm. It is only possible to
% reach a point pi away on a perfect shell. Thus the gauge is limited to an
% arc smaller than a semi circle with a radius smaller than the shell
% radius. For gauges of length lgx = 4*sqrt(rt), this is satisfied for r/t
% > (16/pi^2) i.e. r/t > 1.62 which is very likely the case for in-service
% shells. The other gauge length lg0 is limited to be smaller than r so
% will also be smaller than a semi circle.

searchLimit = find(thetv > thetv(j) & thetv < thetv(j) + pi, 1, 'last'); % < and not <= is used to allow the segment which may cross over into >pi range to be considered.

end

function [snapTwoPos, interPairsMatch, intersects, snapOneX, snapOneY, gaugePointX, gaugePointY, gaugeEndX, gaugeEndY, snapTwoPosX, snapTwoPosY, searchRange] = computeShellGaugeCorrespondences(xv, yv, x1, y1, xcent, ycent, lgxtheta, rnom, searchLimit, snapOneInd, k, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise searchRange, the indexes of shellPoints to check from snapOneInd until the end 
% of the shellPoints.
% Initialise snapOne, the coordinates of the point about which the gauge will rotate.
% Compute snapOne2EndLength, the distance between snapOne and gaugeEnd.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Finding snapOne coordinates
snapOneX = xv(snapOneInd); snapOneY = yv(snapOneInd);
r = rnom(k);

%%% CORRESPONDENCES ALONG LINE SEGMENTS ON WALL %%%
searchRange = snapOneInd:searchLimit; % Using this search range for end intersection 

% Identifying line segments where the end of the gauge can land
[gaugeEndX, gaugeEndY] = gaugePoints(x1, y1, xcent, ycent, lgxtheta, r); % End of gauge coordinates
snapOne2EndLength = sqrt((snapOneX-gaugeEndX)^2 + (snapOneY-gaugeEndY)^2); % Distance between end of gauge and SnapOne

adjacentPtLength = sqrt((gaugeEndX - snapOneX)^2 + (gaugeEndY - snapOneY)^2);
dist2AllPts = sqrt((snapOneX-xv(searchRange)).^2 +(snapOneY-yv(searchRange)).^2); % Distance to all points on the circumference from SnapOne

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find interPairs, the line segments along the shell where the end of the gauge could 
% possibly land. These are stored as pairs of indexes.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing line segments where one end is inside circle of radius
% snapOne2EndLength and origin snapOne and the other end is outside i.e.
% single intersection point.

mlPair = find(dist2AllPts(1:end-1) >= snapOne2EndLength & dist2AllPts(2:end) <= snapOne2EndLength)'; % More than r, then less than r
lmPair = find(dist2AllPts(1:end-1) <= snapOne2EndLength & dist2AllPts(2:end) >= snapOne2EndLength)'; % Less than r, then more than r
interPairs = [mlPair mlPair+1; lmPair lmPair+1;]+snapOneInd-1; % Line segments

% -------------------------------------------------------------------- %
%%% Trying to catch edge case where the ends of the line segment are both
%%% outside the above mentioned circle i.e. further than snapOne2EndLength
%%% away from snapOne but there is still an intersection (1 tangential
%%% intersection or two intersections between line and circle.)

% Possible line segments
interPairExtra = [searchRange(1:end-1); searchRange(2:end)];

% Removing line segments where both ends are inside the mentioned circle
sel = dist2AllPts(1:end-1) < snapOne2EndLength & dist2AllPts(2:end) < snapOne2EndLength;
interPairExtra(:,sel) = [];

% Start and end coordinates of remaining line segments
xvStart = xv(interPairExtra(1,:)); yvStart = yv(interPairExtra(1,:));
xvStartp1 = xv(interPairExtra(2,:)); yvStartp1 = yv(interPairExtra(2,:));

% Normalised vector along line segment
lineVector = [xvStartp1 - xvStart;yvStartp1 - yvStart;];
lineVectorLen = vecnorm([xvStartp1 - xvStart;yvStartp1 - yvStart;],2,1);
lineVector = lineVector./lineVectorLen;
a = lineVector(1,:);
b = lineVector(2,:);

% Computing coordinate of closest point on line segment (dot product of
% vector from centre of circle to line segment and vector along line
% segment must be zero)
lam = -1*(a.*(xvStart-snapOneX) + b.*(yvStart-snapOneY))./(a.^2+b.^2);
closeVect = [(xvStart-snapOneX)+a.*lam; (yvStart-snapOneY)+b.*lam];

% Computing distance to closest point on line segment from centre of circle
closeDist = sqrt(closeVect(1,:).^2 + closeVect(2,:).^2);

% Finding points where the closest point to the centre of the circle is
% also on the line segment. Removing if it is not on line segment.
sel = lam >0 & lam <= lineVectorLen;
interPairExtra = interPairExtra(:,sel)';
closeDist = closeDist(sel);

% Finding points where closest point on the line segment is closer than
% snapOne2EndLength i.e. there is an intersection possible
sel = closeDist <= snapOne2EndLength+64*eps(max(1,snapOne2EndLength));
interPairExtra = interPairExtra(sel,:);
% -------------------------------------------------------------------- %

% Combining both solutions and removing any doubles if any
interPairs = [interPairs; interPairExtra];
interPairs = unique(interPairs, 'rows');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for pair in interPairs do
%     Solve analytically for the intersection between the line segment and a circle of radius 
%     snapOne2EndLength and origin snapOne and store in array intersects.
% end for  
% Append intersects to array snapTwoPosPoints, and gaugeEnd to gaugePoints for each 
% intersection of intersects found if any.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Solving for intersections of gauge with line segments
intersects = zeros(2*size(interPairs,1),2);
interPairsMatch = zeros(2*size(interPairs,1),2); % One-to-one correspondence with intersections
intersectionCount = 0;
for pairInd = 1:size(interPairs,1)
    [posSols] = lineCircleIntersect(interPairs(pairInd,1), snapOneX, snapOneY, adjacentPtLength, xv, yv);
    if ~isempty(posSols)
        target = intersectionCount+(1:size(posSols,1));
        intersects(target,:) = posSols;
        interPairsMatch(target,:) = repmat(interPairs(pairInd,:),size(posSols,1),1);
        intersectionCount = target(end);
    end
end
intersects = intersects(1:intersectionCount,:);
interPairsMatch = interPairsMatch(1:intersectionCount,:);


%%% CORRESPONDENCES WITH VERTICIES ON WALL %%%

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find snapTwoPos, the indexes of points of shellPoints with an index greater than or equal 
% to snapOneInd within snapOne2EndLength of snapOne where the gauge could rotate to.
% Extract coordinates of points along shell wall, snapTwoPosPoints, with index snapTwoPos.
% Compute points, gaugePoints, along the gauge with the same distance from snapOne as 
% snapTwoPosPoints.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Finding relevant points
searchRange = searchRange(2:end);
dist2AllPts = dist2AllPts(2:end);
distanceCheck = dist2AllPts < snapOne2EndLength;
snapTwoPos = searchRange(distanceCheck);

% snapTwo coordinates
snapTwoPosX = xv(snapTwoPos); snapTwoPosY = yv(snapTwoPos);

% Distance and angle enclosed from snapOne to snapTwo
adjacentPtLength = sqrt((snapTwoPosX - snapOneX).^2 + (snapTwoPosY - snapOneY).^2);
beta = 2*asin(max(-1,min(1,adjacentPtLength/(2*r))));

% Defining points along gauge
[gaugePointX, gaugePointY] = gaugePoints(x1, y1, xcent, ycent, beta+atan2(snapOneY-ycent, snapOneX-xcent)-atan2(y1-ycent, x1-xcent), r);


% Adding end of gauge correspondences
gaugePointX = [gaugePointX gaugeEndX*ones(1,size(intersects,1))]; gaugePointY = [gaugePointY gaugeEndY*ones(1,size(intersects,1))];
snapTwoPosX = [snapTwoPosX intersects(:,1)']; snapTwoPosY = [snapTwoPosY intersects(:,2)'];


if plotting
    figure('Units', 'pixels','Position', [1920+680 458 900 500])
    tiledlayout(1,4)
    nexttile
    hold on
    plot(xv, yv, '-xk')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    plot(snapTwoPosX, snapTwoPosY, 'bo')
    % plot(xshell, yshell, 'k','Marker','x', 'MarkerFaceColor','r')
    plot(gaugePointX, gaugePointY, 'mx')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Considered snapping locations', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

end

function [x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, tolMeasureRange, gaugeSnapInd, snapOneInd, snapOneIndOld, endReached] = applyGaugeRotations(gaugePointX, gaugePointY, snapTwoPosX, snapTwoPosY, snapOneInd, snapOneX, snapOneY, snapTwoPos, x1, y1, xcent, ycent, gaugeEndX, gaugeEndY, xv, yv, interPairs, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the angle enclosed by vectors from snapOne to gaugePoints and snapOne to 
% snapTwoPosPoints. 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing angle using cross product
av = [gaugePointX-snapOneX;gaugePointY-snapOneY];
bv = [snapTwoPosX-snapOneX;snapTwoPosY-snapOneY];
crossValue = av(1,:).*bv(2,:)-av(2,:).*bv(1,:);
dotValue = sum(av.*bv,1);
snapRotTwoAng = atan2(crossValue,dotValue);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine the appropriate angle to use based on its sign and whether the inside or outside 
% of the shell is being assessed.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% In the extremely unlikely case there are two points with the exact same
% rotation angle, the one with the smaller index would be picked.

% This algorithm is defined for outside gauge placement.
[~,gaugeSnapInd] = min(snapRotTwoAng);
snapRotTwoAng = snapRotTwoAng(gaugeSnapInd);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Rotate gaugeStart and gaugeEnd about snapOne.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Rotation matrix to facilitate second snap
Rotm = [cos(snapRotTwoAng) -sin(snapRotTwoAng); sin(snapRotTwoAng) cos(snapRotTwoAng)];

% Rotating rod centre
rotated = Rotm*[xcent-snapOneX; ycent-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
xcent = rotated(1,:); ycent = rotated(2,:); % Extracting rotated coordinates

% Rotating rod start
rotated = Rotm*[x1-snapOneX; y1-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
x1 = rotated(1,:); y1 = rotated(2,:); % Extracting rotated coordinates

% Rotating gauge point - Just for plotting
rotated = Rotm*[gaugePointX-snapOneX; gaugePointY-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
gaugePointX = rotated(1,:); gaugePointY = rotated(2,:); % Extracting rotated coordinates

% Rotating gauge end - Just for plotting
rotated = Rotm*[gaugeEndX-snapOneX; gaugeEndY-snapOneY]+[snapOneX; snapOneY]; % Performing rotation
gaugeEndX = rotated(1,:); gaugeEndY = rotated(2,:); % Extracting rotated coordinates

if plotting
    nexttile
    hold on
    plot(xv, yv, '-kx')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    plot(gaugePointX, gaugePointY, 'mx')
    plot(gaugePointX(gaugeSnapInd), gaugePointY(gaugeSnapInd), 'bo')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Snap', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine snapTwo, the second point of contact of the gauge with the shell wall.
% Determine tolMeasureRange, the indexes of points where tolerance measurements are 
% feasible.
% Determine the boolean endReached whether the end of the shell or gauge has been reached.
% if not endReached then
%     Update snapOneInd with index of snapTwo.
% end if
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Next rotation
if gaugeSnapInd <= length(snapTwoPos) % End of gauge not reached
    tolMeasureRange = snapOneInd+1:(snapOneInd + gaugeSnapInd)-1;
    snapOneIndOld = snapOneInd;
    snapOneInd = snapTwoPos(gaugeSnapInd);
    endReached = false;

% End of gauge reached
else
    tolMeasureRange = snapOneInd+1:interPairs(gaugeSnapInd-length(snapTwoPos),1);
    snapOneIndOld = snapOneInd;
    endReached = true;
end

end

function [measureTol] = selfIntersectionCheck(x1, y1, xcent, ycent, snapOneX, snapOneY, xv, yv, rnom, j, k, snapOneIndOld, lgxtheta, plotting)
numCheckPts = 100;
r = rnom(k);

% Compupting test gauge
snapOneTheta = atan2(snapOneY-ycent, snapOneX-xcent); % Position of snapOne relative to gaugeCentre
startTheta = atan2(y1-ycent, x1-xcent); % Position of start of gauge relaltive to gaugeCentre

snapOneThetaRelStart = snapOneTheta-startTheta; % Angular difference between snapOne and gaugeStart
if snapOneThetaRelStart < 0; snapOneThetaRelStart = snapOneThetaRelStart+2*pi; end % Correcting so numbers positive

% Region from gaugeStart to original snapOne
if j == snapOneIndOld % Making sure original snapOne is not same as gaugeStart
    firstRegionCheckPointsTheta = [];
else
    firstRegionCheckPointsTheta = linspace(0,snapOneThetaRelStart,numCheckPts); % Points between gaugeStart and snapOne
end

secondRegionCheckPointsTheta = linspace(snapOneThetaRelStart,lgxtheta,numCheckPts); % Points between snapOne and gaugeEnd
checkPointsTheta = [firstRegionCheckPointsTheta(2:end-1) secondRegionCheckPointsTheta(2:end-1)]; % Removing first and last values of each region to avoid end points on shell wall being considered
[gaugeCheckX, gaugeCheckY] = gaugePoints(x1, y1, xcent, ycent, checkPointsTheta, r); % Cartesian coordinates of points to consider

% Checking for intersections
[in,on] = inpolygon(gaugeCheckX,gaugeCheckY,xv,yv);

sel = ~in|on;
if any(~sel)
    % Intersection found
    measureTol = false;
else
    % No intersection found
    measureTol = true;
end

if plotting
    nexttile;
    hold on
    plot(xv, yv, '-kx')
    % plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xcent, ycent, 'ko')
    % plot(xc, yc, 'g','Marker','x', 'MarkerFaceColor','r')
    plot(gaugeCheckX(~sel), gaugeCheckY(~sel), 'ro')
    plot(gaugeCheckX(sel), gaugeCheckY(sel), 'go')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Intersection check', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

function [U0Tol] = measureTolerance(x1, y1, xcent, ycent, snapOneX, snapOneY, gaugeEndX, gaugeEndY, gaugeSnapInd, tolMeasureRange, xv, yv, interPairs, rnom, k, lgx, lgxtheta, snapTwoPos, snapOneIndOld, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute normalised vectors n, from gaugeCentre to shellPoints.
% Compute possible points on the gauge which are rnom from gaugeCentre along vectors n 
% in both opposing directions.
% Determine whether the computed points truly lie on the gauge using gaugeStart and θl and 
% remove incorrect points.
% Determine the distance between shellPoints and points on the gauge.
% Normalise distances by the gaugeLength to yield the toleranceValues.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

r = rnom(k);
% Identifying feasible tolerance measurement locations (must be within
% lgxtheta).
% For all points identified in range of snapOne and snapTwo, find points on
% the gauge which have a perpendicular line between the gauge and the shell

%%%%%%%%%%%%%%%%%%%% GAUGE FORMULATION CODE %%%%%%%%%%%%%%%%%%%%
% Normalised vector from shell to centroid
perVec = [xv(tolMeasureRange)-xcent; yv(tolMeasureRange)-ycent];
perVec = perVec./sqrt(perVec(1,:).^2+perVec(2,:).^2);

% Coordinates along arc of gauge on both halves of circle where the
% perpendicular from the gauge could originate from
xc = [xcent + r*perVec(1,:); xcent - r*perVec(1,:)];
yc = [ycent + r*perVec(2,:); ycent - r*perVec(2,:)];
tolMeasureRange = [tolMeasureRange; tolMeasureRange];

% Computing angle swept from start of gauge to points on gauge.
gaugeAngle = atan2(yc-ycent, xc-xcent); % Angular position of points on gauge relative to gaugeCentre
gaugeStartAngle = atan2(y1-ycent, x1-xcent); % Angular position of gaugeStart relative to gaugeCentre
arcAngle = mod(gaugeAngle-gaugeStartAngle,2*pi); % Clockwise angle swept from gaugeStart to points on gauge

% Identifying feasible coordinates of points on gauge based on whether the
% perpendicular lands within the range of the gauge.
sel = arcAngle <= lgxtheta;
xc = xc(sel); yc = yc(sel);
tolMeasureRange = tolMeasureRange(sel);

% %%%%%%%%%%%%%%%%%%%% GAUGE FORMULATION CODE %%%%%%%%%%%%%%%%%%%%

% Extracting dimple tolerance - There may be points which satisfy all
% conditions but do not have a direct line of sight from the gauge to the
% point on the shell wall due to the shell wall itself obstructing the line
% of sight. Tolerance measures are still taken in these cases.

radialDist = sqrt((xc-xv(tolMeasureRange)').^2 + (yc-yv(tolMeasureRange)').^2);
U0Tol = radialDist/lgx;

% Plotting
if plotting
    nexttile
    hold on
    plot(xv, yv, '-kx')
    plot(snapOneX, snapOneY, 'or')
    plot(x1, y1, 'ok', 'MarkerFaceColor','k')
    plot(xc, yc, 'ok')
    plot([xv(tolMeasureRange); xc'], [yv(tolMeasureRange); yc'], 'g')

    % Plotting appropriate points along gauge
    if gaugeSnapInd > length(snapTwoPos)
        alongShellIndsPlotting = snapOneIndOld:interPairs(gaugeSnapInd-length(snapTwoPos),1);
    else
        alongShellIndsPlotting = snapOneIndOld:(snapOneIndOld + gaugeSnapInd);
    end

    perVec = [[xv(snapOneIndOld) xv(alongShellIndsPlotting) gaugeEndX]-xcent; [yv(snapOneIndOld) yv(alongShellIndsPlotting) gaugeEndY]-ycent]; % Perpendicular to rod vector
    perVec = perVec./sqrt(perVec(1,:).^2+perVec(2,:).^2); % Normalised Perpendicular to rod vector

    % Coordinates along arc
    gaugeCheckX = [xcent + r*perVec(1,:); xcent - r*perVec(1,:)];
    gaugeCheckY = [ycent + r*perVec(2,:); ycent - r*perVec(2,:)];

    gaugeAngle = atan2(gaugeCheckY-ycent, gaugeCheckX-xcent); % Angular position of points on gauge relative to gaugeCentre
    gaugeStartAngle = atan2(y1-ycent, x1-xcent); % Angular position of gaugeStart relative to gaugeCentre
    arcAngle = mod(gaugeAngle-gaugeStartAngle,2*pi); % Clockwise angle swept from gaugeStart to points on gauge

    % Identifying places where the perpendiculars do not land on the gauge
    sel = arcAngle <= lgxtheta;
    gaugeCheckX = gaugeCheckX(sel); gaugeCheckY = gaugeCheckY(sel);

    plot(gaugeCheckX, gaugeCheckY, 'bx')
    plot(xcent, ycent, 'ko')
    axis equal tight
    grid on
    xlabel('x [mm]', 'Interpreter','latex')
    ylabel('y [mm]', 'Interpreter','latex')
    title('Tolerance extraction', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end

end

function [gaugePointX, gaugePointY] = gaugePoints(x1, y1, xcent, ycent, beta, r)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the angular position corresponding to beta.
% Convert the angular position corresponding to beta into Cartesian coordinates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Defining points along gauge
startAng = atan2(y1-ycent, x1-xcent); % Angular position of (x1, y1) relative to gaugeCentre
endAng = startAng + beta; % Angular position of gaugePoints relative to gaugeCentre

% Coordinates of gaugePoints
gaugePointX = r*cos(endAng) + xcent;
gaugePointY = r*sin(endAng) + ycent;
end

% Solving for the intersection of shell wall and end of gauge
function [sols] = lineCircleIntersect(startInd, snapOneX, snapOneY, adjacentPtLength, xv, yv)

sols = lineCircleSegmentIntersections( ...
    [xv(startInd),yv(startInd)],[xv(startInd+1),yv(startInd+1)], ...
    [snapOneX,snapOneY],adjacentPtLength);

end

% Best fit circle to data
function [xv, yv, thetv] = bestFit(xv, yv)
[centre,~] = bestFitCircle2D(xv,yv);

% Reiorientating geometry
xv = xv-centre(1);
yv = yv-centre(2);

% Compputing circumferential coordinate
thetv = unwrap(atan2(yv,xv));

% Ensuring theta starts at zero contiuously and increases
thetv = thetv-thetv(1);
end
