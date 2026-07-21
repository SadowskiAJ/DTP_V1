% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Code measures the meridional dimple tolerance

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
startInd = 1;
for j = 1:length(STRAKES)
    S = STRAKES(j);
    inds = startInd:N2Cum(j); % Indicies of a particular strake
    % rnom(inds) = req(z(inds), S); % Less accurate
    rnom(inds) = req2(z(inds), S, z(inds(1)), z(inds(end)));
    tnom(inds) = THICKS(S)*ones(size(inds));
    startInd = N2Cum(j)+1;
end

rnom = rnom - tnom/2; % Correcting to obtain nominal shell midsurface
R_MESH = R_MESH - tnom/2; % Converting reconstructed outer surface to actual midsurface

%% Entire surface evaluation
U0x = nan(size(R_MESH)); % Unassessed positions remain explicitly undefined
side = 'outside'; % Choose whether the gauge is to be placed inside or outside the shell
% side = 'inside';
plotting = false; % Whether plots are generated or not. Turn off for looped runs!!!

parfor k = 1:size(R_MESH,2)
    fprintf('Examining k = %g of %g \n', k, size(R_MESH,2))
    rv = R_MESH(:,k)'; zv = Z_MESH(:,k)';
    U0column = nan(size(R_MESH,1),1);
    for j = 1:size(R_MESH,1)-1
        
        % Tolerance measurement for a given start point as described in
        % pseudocode
        [lgx, r1, z1, r2, z2, snapOneInd] = defineMeridionalGauge(rv, zv, rnom, tnom, j, plotting);
        while true
            [snapTwoPos, interPairs, intersects, s1r, s1z, gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, searchRange] = computeShellGaugeCorrespondences(rv, zv, r1, z1, r2, z2, snapOneInd, plotting);
            [r1, z1, r2, z2, tolMeasureRange, snapTwoR, snapTwoZ, snapOneInd, endReached] = applyGaugeRotations(gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, s1r, s1z, snapTwoPos, searchRange, intersects, r1, z1, r2, z2, rv, zv, interPairs, side, plotting);
            [U0Tol] = measureTolerance(r1, z1, r2, z2, rv, zv, s1r, s1z, snapTwoR, snapTwoZ, intersects, tolMeasureRange, lgx, plotting);
            candidate = abs(U0Tol(:));
            current = U0column(tolMeasureRange);
            replace = isfinite(candidate) & ...
                (~isfinite(current) | candidate > current);
            current(replace) = candidate(replace);
            U0column(tolMeasureRange) = current;
            if endReached
                break
            end
        end

    end
    U0x(:,k) = U0column;
end

%% Plotting tolerance variation of U0x 
resultFigure = plotToleranceHeatmap(THET_MESH,Z_MESH,R_MESH,U0x, ...
    'Meridional dimple tolerance $U_{0x}$','$U_{0x}$ [-]');
saveToleranceAssessmentFigure(resultFigure,inputDir, ...
    ['Tolerance_MeridionalDimple_U0x_',tower])

%% Single gauge placement evaluation - Uncomment to see tolerance evaluated at a single gauge start position
% U0xSingle = zeros(size(R_MESH,1),1);
% side = 'outside'; % Choose whether the gauge is to be placed inside or outside the shell
% % side = 'inside';
% plotting = true; % Whether plots are generated or not
% 
% k = 1; % Circumferential index k
% fprintf('Examining k = %g\n', k)
% rv = R_MESH(:,k)'; zv = Z_MESH(:,k)';
% j = 1; % Vertical index j
% 
% 
% [lgx, r1, z1, r2, z2, snapOneInd] = defineMeridionalGauge(rv, zv, rnom, tnom, j, plotting);
% 
% while true
%     [snapTwoPos, interPairs, intersects, s1r, s1z, gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, searchRange] = computeShellGaugeCorrespondences(rv, zv, r1, z1, r2, z2, snapOneInd, plotting);
%     [r1, z1, r2, z2, tolMeasureRange, snapTwoR, snapTwoZ, snapOneInd, endReached] = applyGaugeRotations(gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, s1r, s1z, snapTwoPos, searchRange, intersects, r1, z1, r2, z2, rv, zv, interPairs, side, plotting);
%     [U0Tol] = measureTolerance(r1, z1, r2, z2, rv, zv, s1r, s1z, snapTwoR, snapTwoZ, intersects, tolMeasureRange, lgx, plotting);
%     U0xSingle(tolMeasureRange) = max([U0xSingle(tolMeasureRange) abs(U0Tol')],[],2);
%     if endReached
%         break
%     end
% end


%%
function [lgx, r1, z1, r2, z2, snapOneInd] = defineMeridionalGauge(rv, zv, rnom, tnom, j, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute lgx using the nominal properties at j.
% Initialise gaugeVector = [0 1]T to represent the orientation of a vertical gauge.
% Extract coordinates of gaugeStart from shellPoints.
% Compute coordinates of gaugeEnd.
% Initialise snapOneInd, to j representing the start of the gauge index.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

lgx = 4*sqrt(rnom(j)*tnom(j)); % Length of gauge
rodVec = [0; 1]; % Initial assumed gauge orientation (straight up)

r1 = rv(j); % R Location of start of gauge
z1 = zv(j); % Z Location of start of gauge

r2 = r1 + lgx*rodVec(1); % R Location of end gauge
z2 = z1 + lgx*rodVec(2); % Z Location of end gauge

% Initialising snapOneInd
snapOneInd = j;

if plotting
    figure;
    hold on
    plot(rv, zv, '-x')
    plot(r1, z1, 'ok', 'MarkerFaceColor','m')
    plot(r2, z2, 'ok', 'MarkerFaceColor','y')
    plot([r1 r2], [z1 z2], 'k')
    xlabel('$\rho$ [mm]', 'Interpreter','latex')
    ylabel('z [mm]', 'Interpreter','latex')
    title('Plotting initial gauge starting position', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end

end

function [snapTwoPos, interPairsMatch, intersects, s1r, s1z, gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, searchRange] = computeShellGaugeCorrespondences(rv, zv, r1, z1, r2, z2, snapOneInd, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise searchRange, the indexes of shellPoints to check from snapOneInd until the end 
% of the shellPoints.
% Initialise snapOne, the coordinates of the point about which the gauge will rotate.
% Compute snapOne2EndLength, the distance between snapOne and gaugeEnd.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

searchRange = snapOneInd:length(zv);
s1r = rv(searchRange(1)); s1z = zv(searchRange(1)); % Coordinates about which gauge will rotate

% Snap range based on distance computed for every new SnapOne location
snapOne2EndLength = sqrt((s1r-r2)^2 + (s1z-z2)^2); % Distance between end of gauge and SnapOne
dist2AllPts = sqrt((s1r-rv(searchRange)).^2 +(s1z-zv(searchRange)).^2); % Distance to all points from SnapOne

%%% CORRESPONDENCES WITH VERTICIES ON WALL %%%

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find snapTwoPos, the indexes of points of shellPoints with an index greater than or equal 
% to snapOneInd within snapOne2EndLength of snapOne where the gauge could rotate to.
% Extract coordinates of points along shell wall, snapTwoPosPoints, with index snapTwoPos.
% Compute points, gaugePoints, along the gauge with the same distance from snapOne as 
% snapTwoPosPoints.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing correspondences along shell wall
distanceCheckRange = find(dist2AllPts < snapOne2EndLength);
distanceCheckRange(distanceCheckRange == 1) = [];
snapTwoPos = searchRange(distanceCheckRange);
snapTwoPosR = rv(snapTwoPos); snapTwoPosZ = zv(snapTwoPos);

% Distance from snapOne to snapTwoPos
snapOne2SnapTwoPosLength = sqrt((snapTwoPosR - s1r).^2 + (snapTwoPosZ - s1z).^2);

% Defining coordinates along gauge
rodVec = [r2-r1; z2-z1];
rodVec = rodVec/norm(rodVec);
gaugePointR = s1r + snapOne2SnapTwoPosLength*rodVec(1);
gaugePointZ = s1z + snapOne2SnapTwoPosLength*rodVec(2);

%%% CORRESPONDENCES ALONG LINE SEGMENTS ON WALL %%%

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find interPairs, the line segments along the shell where the end of the gauge could 
% possibly land. These are stored as pairs of indexes.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing line segments where one end is inside circle of radius
% snapOne2EndLength and origin snapOne and the other end is outside i.e.
% single intersection point.
mlPair = find(dist2AllPts(1:end-1) >= snapOne2EndLength & dist2AllPts(2:end) <= snapOne2EndLength)'; % More than snapOne2EndLength, then less than snapOne2EndLength
lmPair = find(dist2AllPts(1:end-1) <= snapOne2EndLength & dist2AllPts(2:end) >= snapOne2EndLength)'; % Less than snapOne2EndLength, then more than snapOne2EndLength
interPairs = [mlPair mlPair+1; lmPair lmPair+1;]+searchRange(1)-1; % Possible line segments

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
rvStart = rv(interPairExtra(1,:)); zvStart = zv(interPairExtra(1,:));
rvStartp1 = rv(interPairExtra(2,:)); zvStartp1 = zv(interPairExtra(2,:));

% Normalised vector along line segment
lineVector = [rvStartp1 - rvStart;zvStartp1 - zvStart;];
lineVectorLen = vecnorm([rvStartp1 - rvStart;zvStartp1 - zvStart;],2,1);
lineVector = lineVector./lineVectorLen;
a = lineVector(1,:);
b = lineVector(2,:);

% Computing coordinate of closest point on line segment (dot product of
% vector from centre of circle to line segment and vector along line
% segment must be zero)
lam = -1*(a.*(rvStart-s1r) + b.*(zvStart-s1z))./(a.^2+b.^2);
closeVect = [(rvStart-s1r)+a.*lam; (zvStart-s1z)+b.*lam];

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

% Computing coordinates of intersection points
intersects = zeros(2*size(interPairs,1),2);
interPairsMatch = zeros(2*size(interPairs,1),2); % One-to-one correspondence with intersections
intersectionCount = 0;
for pairInd = 1:size(interPairs,1)
    [posSols] = lineCircleIntersect(interPairs(pairInd,1), s1r, s1z, snapOne2EndLength, rv, zv);
    if ~isempty(posSols)
        target = intersectionCount+(1:size(posSols,1));
        intersects(target,:) = posSols;
        interPairsMatch(target,:) = repmat(interPairs(pairInd,:),size(posSols,1),1);
        intersectionCount = target(end);
    end
end
intersects = intersects(1:intersectionCount,:);
interPairsMatch = interPairsMatch(1:intersectionCount,:);

% Adding end of gauge correspondences to relevant vectors
if ~isempty(intersects)
    gaugePointR = [gaugePointR r2*ones(1,size(intersects,1))]; gaugePointZ = [gaugePointZ z2*ones(1,size(intersects,1))];
    snapTwoPosR = [snapTwoPosR intersects(:,1)']; snapTwoPosZ = [snapTwoPosZ intersects(:,2)'];
end

if plotting
    figure
    hold on
    plot(rv, zv, '-xk')
    plot(s1r, s1z, 'or')
    plot(r1, z1, 'ok', 'MarkerFaceColor','m')
    plot(r2, z2, 'ok', 'MarkerFaceColor','y')
    plot(snapTwoPosR, snapTwoPosZ, 'bo')
    plot(gaugePointR, gaugePointZ, 'mx')
    if ~isempty(intersects)
        scatter(intersects(:,1), intersects(:,2), 'rx')
    end
    % axis equal tight
    xlabel('$\rho$ [mm]', 'Interpreter','latex')
    ylabel('z [mm]', 'Interpreter','latex')
    grid on
    title('Considered snapping locations', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

function [r1, z1, r2, z2, tolMeasureRange, snapTwoR, snapTwoZ, snapOneInd, endReached] = applyGaugeRotations(gaugePointR, gaugePointZ, snapTwoPosR, snapTwoPosZ, s1r, s1z, snapTwoPos, searchRange, intersects, r1, z1, r2, z2, rv, zv, interPairs, side, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute the angle enclosed by vectors from snapOne to gaugePoints and snapOne to 
% snapTwoPosPoints. 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Computing angle using cross product
av = [gaugePointR-s1r;gaugePointZ-s1z];
bv = [snapTwoPosR-s1r;snapTwoPosZ-s1z];
crossValue = av(1,:).*bv(2,:)-av(2,:).*bv(1,:);
dotValue = sum(av.*bv,1);
snapRotTwoAng = atan2(crossValue,dotValue);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine the appropriate angle to use based on its sign and whether the inside or outside 
% of the shell is being assessed.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

if strcmp(side, 'inside')
    [~, gaugeSnapInd] = max(snapRotTwoAng);
    snapRotTwoAng = snapRotTwoAng(gaugeSnapInd);
elseif strcmp(side, 'outside')
    [~, gaugeSnapInd] = min(snapRotTwoAng);
    snapRotTwoAng = snapRotTwoAng(gaugeSnapInd);
end 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Rotate gaugeStart and gaugeEnd about snapOne.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Rotation matrix to facilitate second snap
Rotm = [cos(snapRotTwoAng) -sin(snapRotTwoAng); sin(snapRotTwoAng) cos(snapRotTwoAng)];

% Rotating gauge start
rotated = Rotm*[r1-s1r; z1-s1z]+[s1r; s1z]; % Performing rotation
r1 = rotated(1,:); z1 = rotated(2,:); % Extracting rotated coordinates

% Rotating gauge end 
rotated = Rotm*[r2-s1r; z2-s1z]+[s1r; s1z]; % Performing rotation
r2 = rotated(1,:); z2 = rotated(2,:); % Extracting rotated coordinates

% Rotating gaugePoints - Just for plotting
rotated = Rotm*[gaugePointR-s1r; gaugePointZ-s1z]+[s1r; s1z]; % Performing rotation
gaugePointR = rotated(1,:); gaugePointZ = rotated(2,:); % Extracting rotated coordinates

if plotting
    figure;
    hold on
    plot(rv, zv, '-x')
    plot(r1, z1, 'ok', 'MarkerFaceColor','m')
    plot(r2, z2, 'ok', 'MarkerFaceColor','y')
    plot([r1 r2], [z1 z2], 'k')
    if ~isempty(intersects)
        plot(intersects(:,1), intersects(:,2),'rx')
    end
    plot(s1r, s1z, 'ro')
    plot(snapTwoPosR(gaugeSnapInd), snapTwoPosZ(gaugeSnapInd), 'go')
    plot(gaugePointR, gaugePointZ, 'mx')
    % axis equal tight
    xlabel('$\rho$ [mm]', 'Interpreter','latex')
    ylabel('z [mm]', 'Interpreter','latex')
    grid on
    title('Performing snap', 'Interpreter','latex')
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

% End of shell reached
if gaugeSnapInd == length(gaugePointR) && size(intersects,1) == 0
    tolMeasureRange = searchRange(2:end-1);
    snapTwoR = rv(end); snapTwoZ = zv(end); % End of shell coordinates
    snapOneInd = nan; % No longer needed
    endReached = true;
    
% Next rotation
elseif gaugeSnapInd <= length(snapTwoPos)
    tolMeasureRange = searchRange(2):(snapTwoPos(gaugeSnapInd)-1);
    snapTwoR = rv(snapTwoPos(gaugeSnapInd)); snapTwoZ = zv(snapTwoPos(gaugeSnapInd)); % Point along shell coordinate
    snapOneInd = snapTwoPos(gaugeSnapInd);
    endReached = false;
    
% End of gauge reached
else
    tolMeasureRange = searchRange(2):interPairs(gaugeSnapInd-length(snapTwoPos),1);
    snapTwoR = r2; snapTwoZ = z2; % End of gauge coordinates
    snapOneInd = nan; % No longer needed
    endReached = true;
end
end

function [U0Tol] = measureTolerance(r1, z1, r2, z2, rv, zv, s1r, s1z, snapTwoR, snapTwoZ, intersects, tolMeasureRange, lgx, plotting)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Determine the perpendicular distance between the shell and the gauge.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

rodVec = [r2-r1; z2-z1];
rodVec = rodVec/norm(rodVec);

% Considered points on shell
ralongShell = rv(tolMeasureRange); zalongShell = zv(tolMeasureRange);
 
% syms a1 b1 z1 r1 c d x y
% S = solve(a1 == r1 + x*c+y*d, b1 == z1 + x*d-y*c);
% Distance along gauge from start to where perpendicular to gauge landing on shell
alongGaugeDist = (ralongShell*rodVec(1) + zalongShell*rodVec(2) - rodVec(1)*r1 - rodVec(2)*z1)/(rodVec(1)^2 + rodVec(2)^2);

% Perpendicular distance between gauge and shell
perpGaugeDist = (ralongShell*rodVec(2) - zalongShell*rodVec(1) - rodVec(2)*r1 + rodVec(1)*z1)/(rodVec(1)^2 + rodVec(2)^2);
perpGaugeDist = abs(perpGaugeDist);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Check points on gauge to see if they truly lie between snapOne and snapTwo and remove 
% incorrect points.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

start2SnapOne = sqrt((r1-s1r).^2 + (z1-s1z).^2); % Distance between start and snapOne
start2SnapTwo = sqrt((r1-snapTwoR).^2 + (z1-snapTwoZ).^2); % Distance between start and snapTwo

% Identifying points where the perpendicular is between snapOne and snapTwo
sel = (alongGaugeDist >= start2SnapOne) & (alongGaugeDist <= start2SnapTwo);
perpGaugeDist(~sel) = nan; % Removing selected points.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Normalise perpendicular distances by the gauge length to yield the toleranceValues.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

U0Tol = perpGaugeDist/lgx;

if plotting
    % Points on gauge where imperfection is measured to. % Not needed for calc.
    % Only for plotting
    rc = r1 + alongGaugeDist*rodVec(1);
    zc = z1 + alongGaugeDist*rodVec(2);
    figure;
    hold on
    plot(rv, zv, '-x')
    plot(r1, z1, 'ok', 'MarkerFaceColor','m')
    plot(r2, z2, 'ok', 'MarkerFaceColor','y')
    plot([r1 r2], [z1 z2], 'k')
    if ~isempty(intersects)
        plot(intersects(:,1), intersects(:,2),'rx')
    end
    plot(s1r, s1z, 'ro')
    % plot(snapTwoPosR(gaugeSnapInd), snapTwoPosZ(gaugeSnapInd), 'go')
    plot(snapTwoR, snapTwoZ, 'go')
    plot([ralongShell(sel); rc(sel)], [zalongShell(sel); zc(sel)], 'k')
    % axis equal tight
    xlabel('$\rho$ [mm]', 'Interpreter','latex')
    ylabel('z [mm]', 'Interpreter','latex')
    grid on
    title('Tolerance extraction', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

% Solving for the intersection of shell wall and end of gauge
function [sols] = lineCircleIntersect(startInd, snapOneX, snapOneY, adjacentPtLength, xv, yv)

sols = lineCircleSegmentIntersections( ...
    [xv(startInd),yv(startInd)],[xv(startInd+1),yv(startInd+1)], ...
    [snapOneX,snapOneY],adjacentPtLength);

end
