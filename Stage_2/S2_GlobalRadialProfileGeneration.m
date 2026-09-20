function [ZMEDIANS, RMEDIANS] = S2_GlobalRadialProfileGeneration(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower, expectedBeadWidth, numCircWindows, subIntervals, plotting, trimFailedEndJoints, profileInteriorRange, retainLongestProfileRun)

% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Function generates global radial profile

if nargin < 11
    plotting = false;
end
if nargin < 12
    trimFailedEndJoints = false;
end
validateattributes(trimFailedEndJoints,{'logical'},{'scalar'})
if nargin < 14
    retainLongestProfileRun = false;
end
validateattributes(retainLongestProfileRun,{'logical'},{'scalar'})
if nargin < 13
    profileInteriorRange = [-Inf Inf];
end
validateattributes(profileInteriorRange,{'numeric'},{'real','vector','numel',2,'nonnan'})
if profileInteriorRange(1) >= profileInteriorRange(2)
    error('S2_GlobalRadialProfileGeneration:InvalidEndRange','End-trimming regions must not overlap. Reduce profileEndTrimDistance.')
end
if subIntervals < 2 || subIntervals ~= floor(subIntervals)
    error('S2_GlobalRadialProfileGeneration:InvalidSubIntervals', ...
        'subIntervals must be an integer of at least two.')
end
if numCircWindows < 1 || numCircWindows ~= floor(numCircWindows)
    error('S2_GlobalRadialProfileGeneration:InvalidCircumferentialWindows', ...
        'numCircWindows must be a positive integer.')
end

% Stitch only consecutive complete windows. Partial processing keeps one
% continuous profile, since radial offsets cannot be transferred across gaps.

%% Loading data
fprintf('Loading data\n')
[X, Y, Z, ~, ~] = postRegDataLoader(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower);

%% Performing best fit cone and setup

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute best fitting cone and apply transformation to centre and orientate cone upright.
% Compute desired global radial profile vertical spacing ps and vertical window height wh.
% Define windowCentreZ, the vertical coordinates of the centres of the vertical windows 
% using wh and nsub.
% Define windowCentreTheta, the centre of the circumferential window using nc.
% Initialise zG and rG to store the global radial profile.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

fprintf('Performing best fit cone correction \n')
[bestFitParams] = bestFitCone(X, Y, Z, minZ, tower, 1e6);
[X, Y, Z, T, R] = bestFitConeTransform(X, Y, Z, minZ, bestFitParams);

% Creating a copy of data to be accessed during loop
ZCOPY = Z;
dataMaxZ = max(Z);
RCOPY = R;
TCOPY = T;

% Relating to hyperparameters 
% Defining window width etc based on size of expected bead width.
desiredPointSpacing = expectedBeadWidth/10; % Spacing between points of global profile
windowWidth = desiredPointSpacing*subIntervals; % Size of a window

% Window setup
% z coordinate of centres of sliding windows which ensure the results can be
% stiched together into one global profile
windowCentreZ = (min(Z) + windowWidth/2):(windowWidth/subIntervals*(subIntervals-1)):(max(Z) -windowWidth/2);

% Adding one extra window centre to list to make sure the end is reached
if numel(windowCentreZ) < 2
    error('S2_GlobalRadialProfileGeneration:InsufficientHeight', ...
        'The selected elevation range is too short for two profile windows.')
end
dz = windowCentreZ(2) - windowCentreZ(1);
windowCentreZ = [windowCentreZ windowCentreZ(end)+dz];

% Circumferential windows
windowEndstheta = linspace(0, 2*pi, numCircWindows+1); % End points
windowCentreTheta = (windowEndstheta(1:end-1)+windowEndstheta(2:end))/2; % Point at centre of circumferential window

% Arrays to store global profile
ZMEDIANS = [];
RMEDIANS = [];
gapWindow = []; % An incomplete window after the profile has started.
skippedWindows = 0;
longestZ = []; longestR = []; % Keep one continuous region; never join across a gap.

%% Looping through vertical sliding windows
tic

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for vertical sliding window in windowCentreZ do
%   Identify pointsVerticalWindow, the points which lie in the vertical window.
%   Initialise arrays θc, zc, rc to store the circumferential, vertical and radial coordinates after 
%   radial correction.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

for j = 1:length(windowCentreZ)


    fprintf('Examining window %g of %g\n', j, length(windowCentreZ))

    % Copying the initial data
    ZCropped = ZCOPY;
    RCropped = RCOPY;
    TCropped = TCOPY;

    % Identifying correct meridional region
    sel = ZCropped < windowCentreZ(j) - windowWidth/2  | ZCropped > windowCentreZ(j) + windowWidth/2;
    TCropped(sel) = []; RCropped(sel) = []; ZCropped(sel) = [];

    % Initialising arrays to store points after radial correction
    rCirc = [];
    zCirc = [];
    tCirc = [];

    %% Looping through circumferential windows

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   for circumferential window in windowCentreTheta do
    % 	    Identify the points of pointsVerticalWindow which lie in the circumferential window.
    %       Subtract median radius of the circumferential window from the radial values inside it.
    %       Append the radially corrected coordinates to θc, zc, rc
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    for k = 1:length(windowCentreTheta)

        % Identifying correct circumferential region
        if k == 1
            sel = TCropped >= windowEndstheta(k) & TCropped <= windowEndstheta(k+1);
        else
            sel = TCropped > windowEndstheta(k) & TCropped <= windowEndstheta(k+1);
        end

        % Correcting radius by median to account for any circumferential variation in radius
        tCirc(end+1:end+length(RCropped(sel))) = TCropped(sel);
        zCirc(end+1:end+length(RCropped(sel))) = ZCropped(sel);
        rCirc(end+1:end+length(RCropped(sel))) = RCropped(sel)-median(RCropped(sel));

        if ~any(sel); fprintf('j = %g, k = %g EMPTY\n', j,k); end
    end

    %% Defining centre of vertical sub intervals

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   Define arrays zs and zm containing the vertical coordinates of the bounds and centres of 
    %   the nsub sub intervals within the vertical window.
    %   Initialise arrays r_m to store median values within a vertical window.
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    subWindowEndsZ = linspace(windowCentreZ(j) - windowWidth/2, windowCentreZ(j) + windowWidth/2, subIntervals+1);
    subWindowCentreZ = mean([subWindowEndsZ(1:end-1); subWindowEndsZ(2:end)],1);

    % Initialising array to store median values of vertical sub intervals
    subWindowMedianR = nan(size(subWindowCentreZ));

    %% Looping through vertical sub intervals

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   for vertical sub interval in vertical window do
    % 	    Identify points within vertical sub interval and store median of their rc values in rm.
    % 	end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    for k = 1:length(subWindowEndsZ)-1
        Z = zCirc;
        R = rCirc;
        T = tCirc;


        % Identifying points in vertical sub interval
        sel = Z < subWindowEndsZ(k) | Z > subWindowEndsZ(k+1);
        T(sel) = []; R(sel) = []; Z(sel) = [];

        subWindowMedianR(k) = median(R);

    end

    % The extra final window may contain subwindows wholly above the cloud.
    % Drop only those; retain the regular spacing and reject gaps within it.
    if j == length(windowCentreZ)
        beyondData = subWindowEndsZ(1:end-1) > dataMaxZ;
        subWindowCentreZ(beyondData) = [];
        subWindowMedianR(beyondData) = [];
    end

    if isempty(subWindowMedianR) || any(isnan(subWindowMedianR))
        if trimFailedEndJoints && retainLongestProfileRun
            if numel(ZMEDIANS) > numel(longestZ)
                longestZ = ZMEDIANS; longestR = RMEDIANS;
            end
            ZMEDIANS = []; RMEDIANS = [];
            skippedWindows = skippedWindows+1;
            continue
        end
        if trimFailedEndJoints
            skippedWindows = skippedWindows+1;
            if subWindowEndsZ(end) <= profileInteriorRange(1)
                % Discard any isolated usable windows in the lower end region.
                ZMEDIANS = []; RMEDIANS = []; gapWindow = [];
            elseif subWindowEndsZ(1) >= profileInteriorRange(2) && isempty(gapWindow)
                % Stop before the upper-end gap, even if usable windows follow.
                break
            end
            if ~isempty(RMEDIANS) && isempty(gapWindow)
                gapWindow = j;
            end
            continue
        end
        error('S2_GlobalRadialProfileGeneration:EmptySubwindow', ...
            ['Vertical profile window %g contains an empty subwindow. ', ...
             'Increase point density or restrict the height range.'], j)
    end
    if ~isempty(gapWindow)
        error('S2_GlobalRadialProfileGeneration:InteriorGap', ...
            ['Incomplete profile window %g near %.3f m is followed by usable window %g. ', ...
             'Only incomplete end windows may be trimmed; interior gaps cannot be stitched.'], ...
            gapWindow,windowCentreZ(gapWindow),j)
    end

    %% Stitching together global profile

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    %   if first vertical window then
    % 	    Append zm to zG
    %       Correct values in rm so that the array starts at zero and then append to rG.
    % 	else
    % 	    Correct values in rm so that the array starts at last value of rG.
    %       Remove the first value of zm and rm and append to zG and rG respectively.
    % 	end if
    % end for
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    if isempty(RMEDIANS) % First retained vertical window
        ZMEDIANS(end+1:end+length(subWindowCentreZ)) = subWindowCentreZ;
        RMEDIANS(end+1:end+length(subWindowMedianR)) = subWindowMedianR - subWindowMedianR(1);
    else % If any other vertical window
        ZMEDIANS(end+1:end+length(subWindowCentreZ)-1) = subWindowCentreZ(2:end);
        RMEDIANS(end+1:end+length(subWindowMedianR)-1) = subWindowMedianR(2:end) - (subWindowMedianR(1)-RMEDIANS(end));
    end

end
toc

if trimFailedEndJoints && retainLongestProfileRun && numel(longestZ) > numel(ZMEDIANS)
    ZMEDIANS = longestZ; RMEDIANS = longestR;
end

%% Removing bins out of range of data which will contain nans
RMEDIANS(ZMEDIANS > maxZ) = [];
ZMEDIANS(ZMEDIANS > maxZ) = [];

if numel(RMEDIANS) < 2 || any(isnan(RMEDIANS))
    error('S2_GlobalRadialProfileGeneration:InvalidProfile', ...
        'No usable continuous profile was generated; review point density and the height range.')
end
if skippedWindows > 0 && retainLongestProfileRun
    warning('S2_GlobalRadialProfileGeneration:PartialProfile', ...
        ['Incomplete coverage: retained the longest continuous profile, %.3f--%.3f m ', ...
         '(cloud extent %.3f--%.3f m). Data outside this profile are excluded from segmentation; ', ...
         'gaps have not been filled.'],ZMEDIANS(1),ZMEDIANS(end),min(ZCOPY),dataMaxZ)
elseif skippedWindows > 0
    fprintf('Stage 2: trimmed end data; retained profile %.3f--%.3f m.\n', ...
        ZMEDIANS(1),ZMEDIANS(end))
end

%% Plotting
if plotting
    figure;
    plot(ZMEDIANS ,RMEDIANS)
    grid on
    xlabel('$z$ [m]', 'interpreter', 'latex')
    ylabel('$\rho$ [m]', 'interpreter', 'latex')
    title('Global radial profile', 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
end
