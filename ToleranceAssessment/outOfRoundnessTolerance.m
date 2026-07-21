% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Code measures the out-of-roundness tolerance

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
radialScale = (R_MESH-tnom/2)./R_MESH;
X_MESH = X_MESH.*radialScale; % Converting reconstructed outer surface to midsurface
Y_MESH = Y_MESH.*radialScale;
R_MESH = R_MESH-tnom/2;

%% Computing Out of roundness

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% for point in (xv, yv) do
%       Compute distance between point and every other point in (xv, yv).
%       Determine maximum distance and store in a 1D array D. This corresponds to the diameter 
%       measured from the point.
% end for
% Extract the smallest and largest diameter in D and store as dmin and dmax respectively.
% Compute the out-of-roundness Ur of the cross-section using dmin, dmax and dnom.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDOCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

D = zeros(size(THET_MESH)); % Creating empty array to store actual diameters

% Working out diameter when starting from every point in the shell
parfor j = 1:size(THET_MESH,1)
    crossSection = [X_MESH(j,:)',Y_MESH(j,:)'];
    D(j,:) = max(pdist2(crossSection,crossSection),[],2)';
end

Ur = (max(D,[],2) - min(D,[],2))./(2*rnom); % Computing out of roundness

% Plotting histogram of Out of Roundness tolerance Ur
histogramFigure = figure('Color','w','Position',[100 100 1000 650]);
hold on
h = histogram(Ur,'normalization','probability', ...
    'FaceColor',[0.20 0.45 0.75],'EdgeColor','none');
p1 = xline(0.007,'Linestyle','--','Linewidth',1.5,'Color',[0.8500 0.3250 0.0980]);
p2 = xline(0.010,'Linestyle','--','Linewidth',1.5,'Color',[0.9290 0.6940 0.1250]);
p3 = xline(0.015,'Linestyle','--','Linewidth',1.5,'Color',[0.4940 0.1840 0.5560]);
xlabel('$U_r$ [-]','interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title('Distribution of out-of-roundness tolerance $U_r$','interpreter','latex')
legend([p1,p2,p3],{'FTQ Class A','FTQ Class B','FTQ Class C'}, ...
    'interpreter','latex','Location','northeast')
xlim([0 max(max(h.BinEdges),0.015)*1.1])
grid on
box on
saveToleranceAssessmentFigure(histogramFigure,inputDir, ...
    ['Tolerance_OutOfRoundness_Histogram_',tower])


% Plotting variation of Out of Roundness tolerance Ur with meridional height z
heightFigure = figure('Color','w','Position',[100 100 1000 650]);
hold on
plot(Ur,z/1000,'k-','Linewidth',1.4)
p1 = xline(0.007,'Linestyle','--','Linewidth',1.5,'Color',[0.8500 0.3250 0.0980]);
p2 = xline(0.010,'Linestyle','--','Linewidth',1.5,'Color',[0.9290 0.6940 0.1250]);
p3 = xline(0.015,'Linestyle','--','Linewidth',1.5,'Color',[0.4940 0.1840 0.5560]);
xlabel('$U_r$ [-]','interpreter','latex')
ylabel('Elevation $z$ [m]','interpreter','latex')
title('Out-of-roundness tolerance $U_r$ with elevation','interpreter','latex')
legend([p1,p2,p3],{'FTQ Class A','FTQ Class B','FTQ Class C'}, ...
    'interpreter','latex','Location','northeast')
xlim([0 max(max(Ur),0.015)*1.1])
grid on
box on
saveToleranceAssessmentFigure(heightFigure,inputDir, ...
    ['Tolerance_OutOfRoundness_WithHeight_',tower])
