function [regParams, bestFitParamsPreReg, iUse, nClouds] = S1_EnhancedRegistration(minZ, maxZ, units, nperc, tower, numIter, cutoffReg, plotting)

% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Function performs enhanced registration

if nargin < 8
    plotting = false;
end
repoRoot = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(repoRoot, ['Input_', tower]);

%% Loading and plotting data

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute number of data clouds, numClouds.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

files = dir(fullfile(inputDir,'*.bin')); % Finds all files which are .bin
nClouds = length(files); % Number of clouds
if nClouds < 2
    error('S1_EnhancedRegistration:InsufficientClouds', ...
        'At least two .bin point clouds are required in %s.', inputDir)
end
COBJS = []; % Array to store S1_CLOUD objects

% Finding minimum shell radius
run(fullfile(inputDir, ['S0_Data_', tower]))
rmin = min([R0_BOTS R0_TOPS])/1e3; % Minimum radius of shell

% Loading data, initial preprocessing, object creation
for j = 1:nClouds % Looping through all files
    if j == 1
        COBJS = S1_CLOUD(fullfile(files(j).folder, files(j).name), j); % Creating a S1_CLOUD object and cloud array
    else
        COBJS(end+1) = S1_CLOUD(fullfile(files(j).folder, files(j).name), j); % Adding to the cloud array
    end

    COBJS(j) = COBJS(j).CLOUDDOWNRND(nperc); % Downsampling clouds
    COBJS(j) = COBJS(j).CLOUDZLIM(minZ,maxZ); % Limiting the range of the cloud
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
end

if plotting
    figure;
    hold on
    for j = 1:nClouds
        scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
    end
    xlabel('$\theta$ [rad]', 'interpreter', 'latex')
    ylabel(['$z$ [', units ,']'], 'interpreter', 'latex')
    zlabel(['$\rho$ [', units ,']'], 'interpreter', 'latex')
    title('Cylindrical coordinates WTST', 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on

    figure;
    hold on
    for j = 1:nClouds
        scatter3(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, 1, 'marker', '.')
    end
    xlabel(['$x$ [', units ,']'], 'interpreter', 'latex')
    ylabel(['$y$ [', units ,']'], 'interpreter', 'latex')
    zlabel(['$z$ [', units ,']'], 'interpreter', 'latex')
    title('Cartesian coordinates WTST', 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end

%% Applying best fit cone correction

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Fit a best-fitting truncated cone to all clouds at once, apply transformation and store 
% parameters in bestFitParamsPreReg.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Isolating desired region
X = vertcat(COBJS.X);
Y = vertcat(COBJS.Y);
Z = vertcat(COBJS.Z);

sel = Z < minZ | Z > maxZ;
X(sel) = []; Y(sel) = []; Z(sel) = [];

% Computing best fit cone parameters for entire shell
[bestFitParamsPreReg] = bestFitCone(X, Y, Z, minZ, tower, 1e6);

% Applying transformation to each scan
for j = 1:nClouds % Looping through all clouds
    [COBJS(j).X, COBJS(j).Y, COBJS(j).Z, COBJS(j).T, COBJS(j).R] = bestFitConeTransform(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, minZ, bestFitParamsPreReg);
    ptCloud = pointCloud([COBJS(j).X(:), COBJS(j).Y(:), COBJS(j).Z(:)]);
    COBJS(j).normals = pcnormals(ptCloud);
end

%% Iterative loop for ICP

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise regParams to store registration parameters.
% for iteration in (1 to numIter + 1) do
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Initialising arrays
registrationSolutions = zeros((nClouds-1)*6, numIter); % Registration constants from each iteration
ptp99 = zeros(numIter+1,1); % Array to store point to plane metric

for j = 1:nClouds
    COBJS(j).XNEW = COBJS(j).X;
    COBJS(j).YNEW = COBJS(j).Y;
    COBJS(j).ZNEW = COBJS(j).Z;
end

for iterCounter=1:numIter+1
    fprintf('Iteration = %g\n', iterCounter)
    for j = 1:nClouds
        COBJS(j).X = reshape(COBJS(j).XNEW, length(COBJS(j).XNEW),1);
        COBJS(j).Y = reshape(COBJS(j).YNEW, length(COBJS(j).XNEW),1);
        COBJS(j).Z = reshape(COBJS(j).ZNEW, length(COBJS(j).XNEW),1);
    end


    %% Finding extent of each cloud

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % 	for cloud in clouds do
    % 		Divide cloud into circumferential bins and identify if there are points in each bin.
    %       Extend bins circumferentially by necessary amount.
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    CLOUDExtent = {}; % Stores the bins the cloud is present in
    cloudConn = []; % Connections between clouds considered i.e. [1 2; 1 3....]
    CONNECTIONS = {}; % All connections between clouds
    CONNECTIONSdist = {};  % Distance between all connections


    thetOvBins = linspace(0, 2*pi, 101); % Points divinding circumference to determine extent of each scan

    for j = 1:nClouds

        thetExtend = cutoffReg/rmin; % Maximum angle around shell that cutoff distance corresponds to
        thetOvBinSize = thetOvBins(2) - thetOvBins(1); % Size of each theta bin when looking for circumferential extent

        numExtend = ceil(thetExtend/thetOvBinSize); % Number bins to extend circumferential extent accounting for cutoff size
        extentCloud = zeros(length(thetOvBins)-1, 1); % Boolean array storing whether a bin has points in it or not

        % Cylindrical coordinates
        T = atan2(COBJS(j).Y, COBJS(j).X);
        R = sqrt(COBJS(j).X.^2 + COBJS(j).Y.^2);
        T(T < 0) = T(T < 0) + 2*pi;

        % Looping through bins and storing whether there are any points in extentCloud
        for k = 1:length(thetOvBins)-1

            if k == 1
                sel = T >= thetOvBins(k) & T <= thetOvBins(k+1);
            else
                sel = T > thetOvBins(k) & T <= thetOvBins(k+1);
            end

            if any(sel)
                extentCloud(k) = 1;
            end

        end

        % Extending extentCloud variable to account for cutoffReg for acceptable pairings in ICP
        extentInds = find(extentCloud); % Index of bins with points in them
        extendedInds = extentInds; % Index of bins after extension by numExtend
        for offset = 1:numExtend
            extendedInds = union(extendedInds,extentInds-offset);
            extendedInds = union(extendedInds,extentInds+offset);
        end

        % Ensuring the bins after extending are within appropriate ranges
        extendedInds(extendedInds <= 0) = extendedInds(extendedInds <= 0) + length(extentCloud);
        extendedInds(extendedInds > length(extentCloud)) = extendedInds(extendedInds > length(extentCloud)) - length(extentCloud);
        extendedInds = unique(extendedInds);

        % Storing the extended extent of each cloud
        CLOUDExtent{end+1} = extendedInds;
    end


    %% Finding all connections to be used in registration

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   for firstCloudIndex in (1 to numClouds) do
    % 		for secondCloudIndex in (firstCloudIndex + 1 to numClouds) do
    % 			for overlappingBins between firstCloudIndex and secondCloudIndex clouds do
    % 				for point in overlappingBin from cloud with secondCloudIndex do
    % 					Find nearest neighbour of point, to points in cloud with firstCloudIndex
    % 				end for
    % 			end for
    %           Remove nearest neighbours with distances more than cutoffReg.
    %           Store correspondences in cell array called connections.
    % 		end for
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    tic
    for j = 1:nClouds
        points = [COBJS(j).X COBJS(j).Y COBJS(j).Z];
        neighbourSearcher = createns(points, 'NSMethod', 'kdtree', 'Distance', 'euclidean');

        for m = j+1:nClouds
            fprintf('1NN for clouds %g and %g \n', j, m)
            cloudConn(end+1,1:2) = [j m]; % Storing cloud pair being considered

            % Finding overlaps between two clouds
            intersectInds = intersect(CLOUDExtent{j}, CLOUDExtent{m});

            % Radial coordinate of cloud m
            T = atan2(COBJS(m).Y, COBJS(m).X);
            T(T < 0) = T(T < 0) + 2*pi;

            allPoints = [];

            % Looping through overlapping bins
            for n = 1:length(intersectInds)

                % Identifying points in bin
                if intersectInds(n) == 1
                    sel = T >= thetOvBins(intersectInds(n)) & T <= thetOvBins(intersectInds(n)+1);
                else
                    sel = T > thetOvBins(intersectInds(n)) & T <= thetOvBins(intersectInds(n)+1);
                end
                searchOverlapInds = find(sel);
                allPoints = [allPoints; searchOverlapInds(:)]; %#ok<AGROW>
            end

            if isempty(allPoints)
                nearestInds = zeros(0,1);
                dists = zeros(0,1);
            else
                queries = [COBJS(m).X(allPoints), COBJS(m).Y(allPoints), COBJS(m).Z(allPoints)];
                [nearestInds, dists] = knnsearch(neighbourSearcher, queries, 'K', 1);
            end
            
            % Removing nearest neighbour connections outside threshold
            sel = dists > cutoffReg;
            nearestInds(sel) = [];
            allPoints(sel) = [];
            dists(sel) = [];

            % Storing feasible connections and distances between points
            CONNECTIONS{end+1} = [nearestInds allPoints];
            CONNECTIONSdist{end+1} = dists;
        end
    end
    toc

    fprintf('Number of connections = %g\n', sum(cellfun(@numel, CONNECTIONSdist)))


    %% Estimating point to plane before registration

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   Compute point-to-plane metric for all correspondences and store 99th percentile.
    %   if iteration equals (numIter+1) then
    % 		break
    %   end if
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % Array to store point to plane metric
    mindists = [];

    % Computing point to plane metric for two clouds
    for j = 1:size(cloudConn,1)
        c1 = cloudConn(j,:);
        c2 = c1(2); c1 = c1(1); % Relavant clouds
        conn = CONNECTIONS{j}; % Connections between clouds
        if ~isempty(conn)
            c1conn = conn(:,1); c2conn = conn(:,2);
            point2Planedists = abs(dot(COBJS(c1).normals(c1conn,:), [COBJS(c1).X(c1conn)-COBJS(c2).X(c2conn) COBJS(c1).Y(c1conn)-COBJS(c2).Y(c2conn) COBJS(c1).Z(c1conn) - COBJS(c2).Z(c2conn)], 2)); % Point to plane
            mindists(end+1:end+length(point2Planedists)) = point2Planedists;
        end
    end

    if isempty(mindists)
        error('S1_EnhancedRegistration:NoCorrespondences', ...
            'No point-cloud correspondences survived the registration cutoff.')
    end
    fprintf('Point to plane distance p99 = %g m\n', prctile(mindists,99))
    ptp99(iterCounter) = prctile(mindists,99);

    % If statement is neccessary as the point to plane metric is computed
    % for iteration numIter on loop counter numIter+1
    if iterCounter == numIter+1
        break
    end

    %% Solving for registration constants

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   Solve ICP minimisation problem with point-to-plane metric using lsqnonlin for all clouds 
    %   simultaneously and obtain transformation parameters.
    %   Apply solved registration transformation to clouds and append to regParams.
    % end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    nopt = 6;

    % lsqnonlin solver settings
    options = optimoptions('lsqnonlin');
    options.FunctionTolerance = 1e-6;
    options.OptimalityTolerance = 1e-6;
    options.StepTolerance = 1e-6;
    options.MaxFunctionEvaluations = 1e4;
    options.MaxIterations = 1e4;
    options.UseParallel = true;

    xm = zeros((nClouds-1)*nopt,1);
    pairData = buildICPPairData(COBJS, cloudConn, CONNECTIONS);
    clear CONNECTIONS CONNECTIONSdist cloudConn
    f = @(b)ICPFunction(b, pairData, nClouds);

    % Appling lsqnonlin solver
    tic
    x = lsqnonlin(f,xm,[],[],options);
    toc

    % Adding transformations to first one
    regResult = [0 0 0 0 0 0 x']';

    % Applying registration
    for j = 1:nClouds % Looping through all clouds
        COBJS(j) = COBJS(j).CLOUDMERGE(regResult); % Applying registration
        COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).XNEW, COBJS(j).YNEW); % Computing cylindrical coordinates
        COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
    end

    % Storing registration result
    registrationSolutions(:,iterCounter) = x;

end

%% Plotting point to plane metric

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Identify iUse, the iteration with smallest 99th percentile of point-to-plane metric.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

if plotting
    figure
    plot(0:length(ptp99)-1, ptp99, '-x')
    xlabel('Iteration [-]', 'Interpreter','latex')
    ylabel('Point to plane distance [m]', 'Interpreter','latex')
    title('Point to plane distance variation with iteration', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end

% Finding iteration with the smallest point to plane metric
[~, iUse] = min(ptp99);
iUse = iUse-1;
fprintf('Most optimum fit is iteration %g\n', iUse)

% Adding zero motion for first cloud;
regParams = [zeros(6,size(registrationSolutions,2)); registrationSolutions];

%% Plotting results
if plotting
% Loading data again
nperc = 0.02;
files = dir(fullfile(inputDir,'*.bin')); % Finds all files which are .bin
COBJS = []; % Array to store cloud objects

% Loading data, initial preprocessing, object creation
for j = 1:nClouds % Looping through all clouds
    if j == 1
        COBJS = S1_CLOUD(fullfile(files(j).folder, files(j).name), j); % Creating a S1_CLOUD object and cloud array
    else
        COBJS(end+1) = S1_CLOUD(fullfile(files(j).folder, files(j).name), j); % Adding to the cloud array
    end

    COBJS(j) = COBJS(j).CLOUDDOWNRND(nperc); % Downsampling clouds
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
end

% Applying original best fit transformation to each scan
for j = 1:nClouds % Looping through all clouds
    [COBJS(j).X, COBJS(j).Y, COBJS(j).Z, COBJS(j).T, COBJS(j).R] = bestFitConeTransform(COBJS(j).X, COBJS(j).Y, COBJS(j).Z, minZ, bestFitParamsPreReg);
end

% Plotting initial geometry
figure
hold on
for j = 1:nClouds % Looping through all clouds
    COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).X, COBJS(j).Y); % Computes cylindrical coordinates
    COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
    scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
end
view(0,0)
grid on
zlim([3.72 3.77])
ylim([29.0 29.1])
xlabel('$\theta$ [rad]', 'interpreter', 'latex')
ylabel('$z$ [m]', 'interpreter', 'latex')
zlabel('$\rho$ [m]', 'interpreter', 'latex')
title(sprintf('Initial geometry @ %g-%g m', gca().YLim(1), gca().YLim(2)), 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')

% Plotting registered geometries
for solNum = 1:numIter
    figure
    hold on
    % Applying registration
    for j = 1:nClouds % Looping through all clouds
        COBJS(j) = COBJS(j).CLOUDMERGEITER(regParams, solNum); % Applying registration
        COBJS(j) = COBJS(j).CLOUDRT(COBJS(j).XNEW, COBJS(j).YNEW); % Computes cylindrical coordinates
        COBJS(j).T(COBJS(j).T<0) = COBJS(j).T(COBJS(j).T<0) + 2*pi;
        scatter3(COBJS(j).T, COBJS(j).Z, COBJS(j).R, 1, 'marker', '.')
    end
    view(0,0)
    grid on
    zlim([3.72 3.77])
    ylim([29.0 29.1])
    xlabel('$\theta$ [rad]', 'interpreter', 'latex')
    ylabel('$z$ [m]', 'interpreter', 'latex')
    zlabel('$\rho$ [m]', 'interpreter', 'latex')
    title(sprintf('Iteration %g @ %g-%g m', solNum, gca().YLim(1), gca().YLim(2)), 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
end
end

end

%% Functions
function pairData = buildICPPairData(clouds, cloudConn, connections)
pairData = cell(size(cloudConn,1),1);
for pairIndex = 1:size(cloudConn,1)
    cloudOne = cloudConn(pairIndex,1);
    cloudTwo = cloudConn(pairIndex,2);
    connection = connections{pairIndex};
    if isempty(connection)
        pairData{pairIndex} = [];
        continue
    end
    indsOne = connection(:,1);
    indsTwo = connection(:,2);
    data.cloudOne = cloudOne;
    data.cloudTwo = cloudTwo;
    data.pointsOne = [clouds(cloudOne).X(indsOne), clouds(cloudOne).Y(indsOne), clouds(cloudOne).Z(indsOne)];
    data.pointsTwo = [clouds(cloudTwo).X(indsTwo), clouds(cloudTwo).Y(indsTwo), clouds(cloudTwo).Z(indsTwo)];
    data.normals = clouds(cloudOne).normals(indsOne,:);
    pairData{pairIndex} = data;
end
end

function residuals = ICPFunction(parameters, pairData, nClouds)
transforms = cell(nClouds,1);
transforms{1} = zeros(6,1);
for cloudIndex = 2:nClouds
    first = (cloudIndex-2)*6+1;
    transforms{cloudIndex} = parameters(first:first+5);
end

counts = cellfun(@(data) pairCount(data), pairData);
residuals = zeros(sum(counts),1);
offset = 0;
for pairIndex = 1:numel(pairData)
    data = pairData{pairIndex};
    if isempty(data)
        continue
    end
    pointsOne = transformPoints(data.pointsOne, transforms{data.cloudOne});
    pointsTwo = transformPoints(data.pointsTwo, transforms{data.cloudTwo});
    count = size(pointsOne,1);
    residuals(offset+(1:count)) = dot(data.normals, pointsOne-pointsTwo, 2);
    offset = offset+count;
end
end

function count = pairCount(data)
if isempty(data)
    count = 0;
else
    count = size(data.pointsOne,1);
end
end

function transformed = transformPoints(points, constants)
Rx = [1 0 0; 0 cos(constants(1)) -sin(constants(1)); 0 sin(constants(1)) cos(constants(1))];
Ry = [cos(constants(2)) 0 sin(constants(2)); 0 1 0; -sin(constants(2)) 0 cos(constants(2))];
Rz = [cos(constants(3)) -sin(constants(3)) 0; sin(constants(3)) cos(constants(3)) 0; 0 0 1];
transformed = points*(Rz*Ry*Rx)' + constants(4:6)';
end

