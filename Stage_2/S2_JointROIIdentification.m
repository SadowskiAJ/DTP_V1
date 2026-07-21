function [jointCoordinates] = S2_JointROIIdentification(ZMEDIANS, RMEDIANS, nJoints, expectedBeadWidth, plotting)

% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Function identifies ROI of joints

if nargin < 5
    plotting = false;
end
if numel(ZMEDIANS) < 2 || numel(ZMEDIANS) ~= numel(RMEDIANS)
    error('S2_JointROIIdentification:InvalidProfile', ...
        'The global radial profile must contain at least two coordinate pairs.')
end

% A gauge is used to work out how much the shell wall deviates from being
% straight. Data is aggregated using moving maximum window. Local maxima
% are found. Local maxima are sorted from largest to smallest and knowledge
% of how many joints there are is used to determine which maxima are joints.
% The largest local maxima are extracted and can be used for joint
% segmentation

%% Compute length of reference line

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Compute length of reference line l using bw.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

lg = expectedBeadWidth*1.5; % Reference line length l

%% Initialising array to store maximum deviations

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Initialise an array dmax to negative infinity to store the maximum deviation for each point 
% along profile.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

peakDeviations = -inf*ones(length(ZMEDIANS), 1);

%% Looping along global profile

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for point of index j in globalRadialProfile do
%   Find the point jn approximately l away from j.
%   Fit a straight line between j and jn.
%   Interpolate the radial coordinates along the line at vertical coordinates defined by zG.
%   Compute the difference between rG and the interpolated radial values along the line. 
%   Find the maximum deviation dmax,j between the line and global radial profile and its 
%   location along the global radial profile jmax.
% 
%   if dmax,j > dmax(jmax) then
% 	    dmax(jmax) = dmax,j
%   end if
% end for
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Spacing between points
dz = ZMEDIANS(2)-ZMEDIANS(1); 

% Number of indexes equal in length to reference line
lgStep = ceil(lg/dz);

for j = 1:length(ZMEDIANS)-lgStep

    % Points on shell
    z1 = ZMEDIANS(j); r1 = RMEDIANS(j);
    z2 = ZMEDIANS(j+lgStep); r2 = RMEDIANS(j+lgStep);

    % Line segment
    zact = ZMEDIANS(j:j+lgStep);
    ract = RMEDIANS(j:j+lgStep);
    rint = interp1([z1 z2],[r1 r2], zact);

    % Difference between line segment and actual values
    rdiff = ract-rint;

    % Finding maximum value in array
    [M,I] = max(rdiff);

    % Storing maximum value in large array if greater than whatever value is
    % already there.
    peakDeviations(j+I-1) = max([peakDeviations(j+I-1) M]);
end


if plotting
    figure
    plot(ZMEDIANS, peakDeviations)
    xlabel('$z$ [m]', 'Interpreter','latex')
    ylabel('Peak deviations [m]', 'Interpreter','latex')
    title('Variation of peak deviation with height', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end

%% Applying moving maximum filter and plotting maxima

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Apply a moving maximum filter to dmax.
% Determine points of local maxima of dmax.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

peakDeviations = movmax(peakDeviations,lgStep);
TF = islocalmax(peakDeviations);

if plotting
    figure
    hold on
    plot(ZMEDIANS, peakDeviations)
    plot(ZMEDIANS(TF), peakDeviations(TF), 'x')
    xlabel('$z$ [m]', 'Interpreter','latex')
    ylabel('Peak deviations [m]', 'Interpreter','latex')
    title('Moving maximum filter applied', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end


%% Identifying joints

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Find jointCoordinates, the vertical coordinate of the nJoints largest local maxima of dmax.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Finding local maxima and sorting into descending order
locMaxVals = peakDeviations(TF);
locMaxZ = ZMEDIANS(TF);
[~,I] = sort(locMaxVals, 'descend');
if numel(I) < nJoints
    error('S2_JointROIIdentification:InsufficientPeaks', ...
        'Requested %g joints but only %g local maxima were identified.', nJoints, numel(I))
end

% Allocating first nJoints as the joint centres
jointCoordinates = sort(locMaxZ(I(1:nJoints)));

if plotting
    figure
    hold on
    plot(ZMEDIANS, peakDeviations)
    plot(locMaxZ(I(1:nJoints)), locMaxVals(I(1:nJoints)), 'x')
    xlabel('$z$ [m]', 'Interpreter','latex')
    ylabel('Peak deviations [m]', 'Interpreter','latex')
    title('Identified joints', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
end


