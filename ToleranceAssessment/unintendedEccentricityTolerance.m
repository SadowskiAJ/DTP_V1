% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Code measures the unintended eccentricity tolerance

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
if any(diff(STRAKES) ~= 1)
    error('unintendedEccentricityTolerance:NoncontiguousStrakes', ...
        'Unintended eccentricity requires adjacent meshed strakes; gaps were detected.')
end
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
R_MESH = R_MESH - tnom/2; % Correcting to obtain real shell midsurface

%% Computing Ue
eaM = zeros(length(STRAKES)-1, length(theta)); % Array to store unintended eccentricity [mm]
etotM = zeros(length(STRAKES)-1, length(theta)); % Array to store total eccentricity [mm]
UeM = zeros(length(STRAKES)-1, length(theta)); % Array to store Unintended eccentricity tolerance [-]
jointElevation = zeros(length(STRAKES)-1,1); % Elevation of each assessed joint [mm]

for j = 1:length(STRAKES)-1
    indBot = N2Cum(j)+1; % Bottom of strake above joint ind
    indTop = N2Cum(j); % Top of strake below joint ind
    jointElevation(j) = mean(z([indTop indBot]));

    eint = rnom(indTop) - rnom(indBot); % Intended eccentricity
    etot = R_MESH(indTop,:) - R_MESH(indBot,:); % Total eccentricity
    ea = etot-eint; % Unintended eccentricity

    eaM(j,1:length(theta)) = ea; % Storing unintended eccentricty
    etotM(j,1:length(theta)) = etot; % Storing Total eccentricty

    tav = mean([tnom(indBot) tnom(indTop)]); % Average thickness at a joint

    Ue = ea/tav; % Ue tolerance
    UeM(j,1:length(theta)) = Ue; % Storing tolerance
    
end

%% Plotting unintended eccentricy tolerance [-]
Uemrs = reshape(UeM,size(UeM,1)*size(UeM,2),1);

LW = 1.1; % Line width
classA = 0.14;
classB = 0.2;
classC = 0.3;
varName = 'U_e';

histogramFigure = figure('Color','w','Position',[100 100 1000 650]);
hold on
h = histogram(abs(Uemrs),'normalization','probability', ...
    'FaceColor',[0.20 0.45 0.75],'EdgeColor','none');
p1 = xline(classA,'color',[0.8500 0.3250 0.0980],'Linestyle','--', 'Linewidth',LW);
p2 = xline(classB,'color',[0.9290 0.6940 0.1250],'Linestyle','--', 'Linewidth',LW);
p3 = xline(classC,'color',[0.4940 0.1840 0.5560],'Linestyle','--', 'Linewidth',LW);
xlabel(sprintf('$|%s|$ [-]',varName),'interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title(sprintf('Distribution of unintended eccentricity tolerance $|%s|$',varName),'interpreter','latex')
legend([p1,p2,p3],{'FTQ Class A','FTQ Class B','FTQ Class C'}, ...
    'interpreter','latex','Location','northeast')
xlim([0 max(max(h.BinEdges),classC)*1.1])
grid on
box on
saveToleranceAssessmentFigure(histogramFigure,inputDir, ...
    ['Tolerance_UnintendedEccentricity_Histogram_',tower])

ecdfFigure = figure('Color','w','Position',[100 100 1000 650]);
hold on
h = cdfplot(abs(Uemrs));
h.LineWidth = LW;
p1 = xline(classA,'color',[0.8500 0.3250 0.0980],'Linestyle','--', 'Linewidth',LW);
p2 = xline(classB,'color',[0.9290 0.6940 0.1250],'Linestyle','--', 'Linewidth',LW);
p3 = xline(classC,'color',[0.4940 0.1840 0.5560],'Linestyle','--', 'Linewidth',LW);
xlabel(sprintf('$|%s|$ [-]',varName),'interpreter','latex')
ylabel('Probability [-]','interpreter','latex')
title(sprintf('ECDF of unintended eccentricity tolerance $|%s|$', varName),'interpreter','latex')
legend([p1,p2,p3],{'FTQ Class A','FTQ Class B','FTQ Class C'}, ...
    'interpreter','latex','Location','southeast')
xlim([0 classC*1.1])
grid on
box on
saveToleranceAssessmentFigure(ecdfFigure,inputDir, ...
    ['Tolerance_UnintendedEccentricity_ECDF_',tower])

jointFigure = figure('Color','w','Position',[100 100 1450 720]);
layout = tiledlayout(jointFigure,1,2,'TileSpacing','compact', ...
    'Padding','compact');
title(layout,'Unintended eccentricity tolerance at circumferential joints', ...
    'Interpreter','latex')

jointAxis = nexttile(layout,1);
imagesc(jointAxis,theta,1:size(UeM,1),abs(UeM))
axis(jointAxis,'xy')
xlim(jointAxis,[0 2*pi])
ylim(jointAxis,[0.5 size(UeM,1)+0.5])
xticks(jointAxis,[0 pi/2 pi 3*pi/2 2*pi])
xticklabels(jointAxis,{'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})
yticks(jointAxis,1:size(UeM,1))
yticklabels(jointAxis,compose('%.1f',jointElevation/1000))
xlabel(jointAxis,'Circumferential coordinate $\theta$ [rad]', ...
    'Interpreter','latex')
ylabel(jointAxis,'Joint elevation $z$ [m]','Interpreter','latex')
title(jointAxis,'Circumferential distribution','Interpreter','latex')
colourBar = colorbar(jointAxis);
colourBar.Label.String = '$|U_e|$ [-]';
colourBar.Label.Interpreter = 'latex';
colormap(jointAxis,turbo)
box(jointAxis,'on')

maximumUe = max(abs(UeM),[],2,'omitnan');
envelopeAxis = nexttile(layout,2);
hold(envelopeAxis,'on')
plot(envelopeAxis,maximumUe,1:size(UeM,1),'o-', ...
    'Color',[0.15 0.15 0.15],'MarkerFaceColor',[0.20 0.45 0.75], ...
    'LineWidth',LW)
p1 = xline(envelopeAxis,classA,'--','Color',[0.8500 0.3250 0.0980], ...
    'LineWidth',LW);
p2 = xline(envelopeAxis,classB,'--','Color',[0.9290 0.6940 0.1250], ...
    'LineWidth',LW);
p3 = xline(envelopeAxis,classC,'--','Color',[0.4940 0.1840 0.5560], ...
    'LineWidth',LW);
ylim(envelopeAxis,[0.5 size(UeM,1)+0.5])
yticks(envelopeAxis,1:size(UeM,1))
yticklabels(envelopeAxis,compose('%.1f',jointElevation/1000))
xlim(envelopeAxis,[0 max([maximumUe; classC])*1.1])
xlabel(envelopeAxis,'Maximum circumferential $|U_e|$ [-]', ...
    'Interpreter','latex')
ylabel(envelopeAxis,'Joint elevation $z$ [m]','Interpreter','latex')
title(envelopeAxis,'Joint envelope','Interpreter','latex')
legend(envelopeAxis,[p1 p2 p3], ...
    {'FTQ Class A','FTQ Class B','FTQ Class C'}, ...
    'Interpreter','latex','Location','southeast')
grid(envelopeAxis,'on')
box(envelopeAxis,'on')

saveToleranceAssessmentFigure(jointFigure,inputDir, ...
    ['Tolerance_UnintendedEccentricity_Surface_',tower])
