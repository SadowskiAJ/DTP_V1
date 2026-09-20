% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function generateProtocolDiagnostics(tower,minZ,maxZ,force,stages)
%GENERATEPROTOCOLDIAGNOSTICS Save standard protocol diagnostics beside inputs.
%
% Both interactive MATLAB FIG files and PNG previews are written to the
% Input_<tower> directory. Existing diagnostics are regenerated only when
% their source MAT file is newer, unless force is true.
% Optional stages selects diagnostics without requiring later-stage caches.

if nargin < 4
    force = false;
end
if nargin < 5
    stages = 1:6;
end
if ~isnumeric(stages) || isempty(stages) || ~all(ismember(stages(:),1:6))
    error('generateProtocolDiagnostics:InvalidStages','stages must contain values from 1 to 6.')
end
if ~(ischar(tower) || (isstring(tower) && isscalar(tower)))
    error('generateProtocolDiagnostics:InvalidTower','tower must be text.')
end
if ~(isscalar(minZ) && isscalar(maxZ) && isfinite(minZ) && isfinite(maxZ) && maxZ > minZ)
    error('generateProtocolDiagnostics:InvalidElevationRange', ...
        'minZ and maxZ must define a finite increasing elevation range.')
end
if ~(isscalar(force) && (islogical(force) || isnumeric(force)))
    error('generateProtocolDiagnostics:InvalidForce','force must be scalar logical.')
end
force = logical(force);
tower = char(tower);

repoRoot = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(repoRoot,['Input_',tower]);
addpath(fullfile(repoRoot,'Stage_1'))
addpath(fullfile(repoRoot,'Stage_3'))
addpath(fullfile(repoRoot,'Stage_4'))
if ~isfolder(inputDir)
    error('generateProtocolDiagnostics:MissingInputDirectory', ...
        'Input directory does not exist: %s',inputDir)
end

s1Cache = fullfile(inputDir,['S1_Registration_',tower,'.mat']);
s2Cache = fullfile(inputDir,['S2_CanSegmentation_',tower,'.mat']);
s5Cache = fullfile(inputDir,['S5_Surface_',tower,'.mat']);
s6Cache = fullfile(inputDir,['S6_Mesh_',tower,'.mat']);
s2PlotCache = s2Cache;
detectionCache = fullfile(inputDir,['S2_JointDetection_',tower,'.mat']);
if isfile(detectionCache) && (~isfile(s2Cache) || dir(detectionCache).datenum >= dir(s2Cache).datenum)
    s2PlotCache = detectionCache;
end
cacheDependencies = {{s1Cache,s2Cache},{s2PlotCache},{s1Cache,s2Cache,s5Cache}, ...
    {s1Cache,s2Cache,s5Cache},{s1Cache,s5Cache},{s6Cache}};
requiredCaches = unique([cacheDependencies{stages}]);
missing = requiredCaches(~cellfun(@isfile,requiredCaches));
if ~isempty(missing)
    error('generateProtocolDiagnostics:MissingCache', ...
        'Required protocol output is missing: %s',missing{1})
end

if ismember(1,stages) && needsRefresh(inputDir,diagnosticStem(1,tower),{s1Cache,s2Cache},force)
    plotRegistrationDiagnostic(inputDir,tower,minZ,maxZ,s1Cache,s2Cache)
end
if ismember(2,stages) && needsRefresh(inputDir,diagnosticStem(2,tower),s2PlotCache,force)
    plotSegmentationDiagnostic(inputDir,tower,s2PlotCache)
end
refreshStage3 = ismember(3,stages) && needsRefresh(inputDir,diagnosticStem(3,tower),s5Cache,force);
refreshStage4 = ismember(4,stages) && needsRefresh(inputDir,diagnosticStem(4,tower),s5Cache,force);
if refreshStage3 || refreshStage4
    stageDiagnostics = loadOrBuildIntermediateDiagnostics( ...
        tower,minZ,maxZ,s1Cache,s2Cache,s5Cache);
    if refreshStage3
        plotOutlierRemovalDiagnostic(inputDir,tower,stageDiagnostics.stage3)
    end
    if refreshStage4
        plotGriddedSurfaceDiagnostic(inputDir,tower,stageDiagnostics.stage4)
    end
end
if ismember(5,stages) && needsRefresh(inputDir,diagnosticStem(5,tower),s5Cache,force)
    plotSurfaceDiagnostic(inputDir,tower,s5Cache)
end
pointCloudStem = sprintf('DTP_Diagnostic_Stage5_PointCloud_%s',tower);
pointFiles = dir(fullfile(inputDir,'*.bin'));
pointSources = [{s1Cache,s5Cache}, ...
    arrayfun(@(file) fullfile(file.folder,file.name),pointFiles,'UniformOutput',false)'];
if ismember(5,stages) && needsRefresh(inputDir,pointCloudStem,pointSources,force)
    plotPointCloudDiagnostic(inputDir,minZ,maxZ,s1Cache,s5Cache,pointCloudStem)
end
if ismember(6,stages) && needsRefresh(inputDir,diagnosticStem(6,tower),s6Cache,force)
    plotMeshDiagnostic(inputDir,tower,s6Cache)
end
end

function tf = needsRefresh(inputDir,stem,sourceFiles,force)
if force
    tf = true;
    return
end
figFile = fullfile(inputDir,[stem,'.fig']);
pngFile = fullfile(inputDir,[stem,'.png']);
if ~isfile(figFile) || ~isfile(pngFile)
    tf = true;
    return
end
if ischar(sourceFiles) || (isstring(sourceFiles) && isscalar(sourceFiles))
    sourceFiles = cellstr(sourceFiles);
end
sourceTimes = cellfun(@(file) dir(file).datenum,sourceFiles);
% Refresh captions and plotting changes as well as numerical data.
codeInfo = dir([mfilename('fullpath'),'.m']);
sourceTimes(end+1) = codeInfo.datenum;
figInfo = dir(figFile);
pngInfo = dir(pngFile);
tf = min(figInfo.datenum,pngInfo.datenum) < max(sourceTimes);
end

function stem = diagnosticStem(stage,tower)
switch stage
    case 1
        descriptor = 'Registration';
    case 2
        descriptor = 'Segmentation';
    case 3
        descriptor = 'OutlierRemoval';
    case 4
        descriptor = 'GriddedSurface';
    case 5
        descriptor = 'Surface';
    case 6
        descriptor = 'Mesh';
    otherwise
        error('generateProtocolDiagnostics:InvalidStage','Unsupported diagnostic stage.')
end
stem = sprintf('DTP_Diagnostic_Stage%g_%s_%s',stage,descriptor,tower);
end

function plotRegistrationDiagnostic(inputDir,tower,minZ,maxZ,s1Cache,s2Cache)
fprintf('Generating Stage 1 registration diagnostic.\n')
s1 = load(s1Cache,'regParams','bestFitParamsPreReg','iUse');
s2 = load(s2Cache,'jointCoordinates');
files = dir(fullfile(inputDir,'*.bin'));
nClouds = numel(files);
if nClouds < 2
    error('generateProtocolDiagnostics:InsufficientClouds', ...
        'At least two BIN clouds are needed for the registration diagnostic.')
end

sliceElevation = chooseSliceElevation(minZ,maxZ,s2.jointCoordinates);
sliceHalfHeight = max(0.05,0.0025*(maxZ-minZ));
sampleFraction = 0.001;
preTheta = cell(nClouds,1); preRadius = cell(nClouds,1); preZ = cell(nClouds,1);
postTheta = cell(nClouds,1); postRadius = cell(nClouds,1); postZ = cell(nClouds,1);
for cloudIndex = 1:nClouds
    fprintf('  Sampling registration cloud %g of %g.\n',cloudIndex,nClouds)
    [x,y,z] = readBinaryPointSample( ...
        fullfile(files(cloudIndex).folder,files(cloudIndex).name),sampleFraction);
    [x,y,z] = bestFitConeTransform( ...
        x,y,z,minZ,s1.bestFitParamsPreReg,false);
    preTheta{cloudIndex} = mod(atan2(y,x),2*pi);
    preRadius{cloudIndex} = hypot(x,y);
    preZ{cloudIndex} = z;
    if s1.iUse > 0
        [x,y,z] = applyRegistrationTransform( ...
            x,y,z,s1.regParams,cloudIndex,s1.iUse);
        postTheta{cloudIndex} = mod(atan2(y,x),2*pi);
        postRadius{cloudIndex} = hypot(x,y);
        postZ{cloudIndex} = z;
    else
        postTheta{cloudIndex} = preTheta{cloudIndex};
        postRadius{cloudIndex} = preRadius{cloudIndex};
        postZ{cloudIndex} = preZ{cloudIndex};
    end
end

figureHandle = figure('Visible','off','Color','w','Position',[100 100 1400 560]);
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact','Padding','compact');
colours = lines(nClouds);
axesHandles = gobjects(2,1);
for panel = 1:2
    axesHandle = nexttile(layout);
    axesHandles(panel) = axesHandle;
    hold(axesHandle,'on')
    if panel == 1
        thetaParts = preTheta; radiusParts = preRadius; zParts = preZ;
        panelTitle = 'Before enhanced registration';
    else
        thetaParts = postTheta; radiusParts = postRadius; zParts = postZ;
        panelTitle = sprintf('After enhanced registration (transformation %g)',s1.iUse);
    end
    selectedRadii = [];
    for cloudIndex = 1:nClouds
        selection = abs(zParts{cloudIndex}-sliceElevation) <= sliceHalfHeight;
        selected = radiusParts{cloudIndex}(selection);
        selectedRadii = [selectedRadii; selected(:)]; %#ok<AGROW>
    end
    if isempty(selectedRadii)
        close(figureHandle)
        error('generateProtocolDiagnostics:EmptyRegistrationSlice', ...
            'No sampled points were found near z = %g.',sliceElevation)
    end
    referenceRadius = median(selectedRadii);
    for cloudIndex = 1:nClouds
        selection = abs(zParts{cloudIndex}-sliceElevation) <= sliceHalfHeight;
        scatter(axesHandle,thetaParts{cloudIndex}(selection), ...
            1000*(radiusParts{cloudIndex}(selection)-referenceRadius),5, ...
            colours(cloudIndex,:),'filled','MarkerFaceAlpha',0.45)
    end
    grid(axesHandle,'on'); box(axesHandle,'on'); xlim(axesHandle,[0 2*pi])
    xlabel(axesHandle,'$\theta$ [rad]','Interpreter','latex')
    ylabel(axesHandle,'Radial deviation from slice median [mm]','Interpreter','latex')
    title(axesHandle,panelTitle,'Interpreter','latex')
    styleAxes(axesHandle)
end
beforeLimits = ylim(axesHandles(1));
afterLimits = ylim(axesHandles(2));
commonLimits = [min(beforeLimits(1),afterLimits(1)), ...
    max(beforeLimits(2),afterLimits(2))];
ylim(axesHandles(1),commonLimits)
ylim(axesHandles(2),commonLimits)
title(layout,sprintf('Registration comparison near z = %.2f m',sliceElevation), ...
    'Color','k','Interpreter','latex')
diagnosticCaption(layout,{'Colours identify scans in this horizontal slice.', ...
    'Closer agreement suggests better alignment here; this is not a whole-cloud connectivity check.'})
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(1,tower))
end

function [x,y,z] = readBinaryPointSample(fileName,fraction,maximumCount)
if nargin < 3
    maximumCount = 50000;
end
fileID = fopen(fileName,'rb','ieee-le');
if fileID < 0
    error('generateProtocolDiagnostics:FileOpenFailed','Unable to open %s.',fileName)
end
cleaner = onCleanup(@() fclose(fileID));
pointCount = fread(fileID,1,'uint64=>double');
if isempty(pointCount) || pointCount < 1 || pointCount ~= floor(pointCount)
    error('generateProtocolDiagnostics:InvalidBinaryHeader', ...
        'Invalid point count in %s.',fileName)
end
sampleCount = min(maximumCount,max(5000,ceil(pointCount*fraction)));
sampleCount = min(sampleCount,pointCount);
sampleIndices = unique(round(linspace(1,pointCount,sampleCount)));
mapping = memmapfile(fileName,'Offset',8,'Format', ...
    {'double',[3 pointCount],'coordinates'},'Writable',false);
coordinates = mapping.Data.coordinates(:,sampleIndices)';
x = coordinates(:,1); y = coordinates(:,2); z = coordinates(:,3);
end

function [x,y,z] = applyRegistrationTransform(x,y,z,parameters,cloudIndex,iterations)
for solution = 1:iterations
    values = parameters(6*cloudIndex-5:6*cloudIndex,solution);
    rotationX = [1 0 0;0 cos(values(1)) -sin(values(1));0 sin(values(1)) cos(values(1))];
    rotationY = [cos(values(2)) 0 sin(values(2));0 1 0;-sin(values(2)) 0 cos(values(2))];
    rotationZ = [cos(values(3)) -sin(values(3)) 0;sin(values(3)) cos(values(3)) 0;0 0 1];
    transformed = rotationZ*rotationY*rotationX*[x(:)';y(:)';z(:)'] + values(4:6);
    x = transformed(1,:)'; y = transformed(2,:)'; z = transformed(3,:)';
end
end

function elevation = chooseSliceElevation(minZ,maxZ,jointCoordinates)
joints = jointCoordinates(:);
joints = sort(joints(isfinite(joints) & joints > minZ & joints < maxZ));
edges = [minZ; joints; maxZ];
midRange = (minZ+maxZ)/2;
interval = find(edges(1:end-1) <= midRange & edges(2:end) >= midRange,1);
if isempty(interval)
    [~,interval] = max(diff(edges));
end
elevation = mean(edges(interval:interval+1));
end

function plotSegmentationDiagnostic(inputDir,tower,s2Cache)
fprintf('Generating Stage 2 segmentation diagnostic.\n')
s2 = load(s2Cache);
z = s2.ZMEDIANS(:);
r = s2.RMEDIANS(:);
figureHandle = figure('Visible','off','Color','w','Position',[100 100 1200 650]);
axesHandle = axes(figureHandle);
profileHandle = plot(axesHandle,z,1000*(r-median(r,'omitnan')), ...
    'k-','LineWidth',1); hold(axesHandle,'on')
joints = s2.jointCoordinates(:);
joints = joints(isfinite(joints) & joints >= min(z) & joints <= max(z));
if isfield(s2,'JOINTBOUNDS')
    bounds = s2.JOINTBOUNDS(:);
else
    bounds = [];
end
bounds = bounds(isfinite(bounds) & bounds >= min(z) & bounds <= max(z));
centreHandle = plot(axesHandle,nan,nan,'-','Color',[0.8 0.1 0.1]);
boundHandle = plot(axesHandle,nan,nan,':','Color',[0.2 0.4 0.8]);
legendHandles = [profileHandle,centreHandle];
legendLabels = {'Global radial profile','Detected joint centres'};
if isfield(s2,'retainedJoints') && any(~s2.retainedJoints) && nnz(s2.retainedJoints) >= 2
    retained = s2.retainedJoints(:);
    allJointZ = s2.jointCoordinates(:);
    runStarts = find(retained & [true;~retained(1:end-1)]);
    runEnds = find(retained & [~retained(2:end);true]);
    limits = ylim(axesHandle);
    endRegions = [[min(z);allJointZ(runEnds)], [allJointZ(runStarts);max(z)]];
    for region = 1:size(endRegions,1)
        x = endRegions(region,:);
        patch(axesHandle,x([1 2 2 1]),limits([1 1 2 2]),[0.8 0.8 0.8], ...
            'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    end
    excludedHandle = plot(axesHandle,nan,nan,'s','Color',[0.6 0.6 0.6]);
    legendHandles(end+1) = excludedHandle;
    legendLabels{end+1} = 'Outside retained joint runs (excluded)';
end
if isfield(s2,'nominalJointZ')
    nominalHandle = plot(axesHandle,nan,nan,'--','Color',[0.3 0.55 0.3]);
    for elevation = s2.nominalJointZ(:)'
        xline(axesHandle,elevation,'--','Color',[0.3 0.55 0.3],'LineWidth',0.8);
    end
    legendHandles(end+1) = nominalHandle;
    legendLabels{end+1} = 'Nominal junctions';
end
for index = 1:numel(joints)
    xline(axesHandle,joints(index),'Color',[0.8 0.1 0.1],'LineWidth',0.8);
end
for index = 1:numel(bounds)
    xline(axesHandle,bounds(index),'Color',[0.2 0.4 0.8], ...
        'LineStyle',':','LineWidth',0.6);
end
grid(axesHandle,'on'); box(axesHandle,'on')
xlabel(axesHandle,'Elevation $z$ [m]','Interpreter','latex')
ylabel(axesHandle,'Radial profile relative to median [mm]','Interpreter','latex')
title(axesHandle,'Stage 2 global radial profile and detected joints','Interpreter','latex')
if isfield(s2,'JOINTBOUNDS')
    diagnosticCaption(axesHandle,{'Red: detected centres. Blue dotted: fitted bounds. Green dashed: nominal junctions.', ...
        'Check that the marked joints follow the features in the measured profile.'})
    legendHandles(end+1) = boundHandle;
    legendLabels{end+1} = 'Joint bounds';
else
    title(axesHandle,'Stage 2 joint detections before validation','Interpreter','latex')
    diagnosticCaption(axesHandle,{'Red: detected centres. Green dashed: nominal junctions.', ...
        'These detections are unvalidated; joint bounds have not been fitted.'})
end
legendHandle = legend(axesHandle,legendHandles,legendLabels,'Location','best', ...
    'Interpreter','latex');
set(legendHandle,'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25])
styleAxes(axesHandle)
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(2,tower))
end

function stageDiagnostics = loadOrBuildIntermediateDiagnostics( ...
    tower,minZ,maxZ,s1Cache,s2Cache,s5Cache)
cacheVariables = whos('-file',s5Cache);
if any(strcmp({cacheVariables.name},'stageDiagnostics'))
    cached = load(s5Cache,'stageDiagnostics');
    if isfield(cached.stageDiagnostics,'stage3') && ...
            isfield(cached.stageDiagnostics.stage3,'thetaBefore') && ...
            isfield(cached.stageDiagnostics,'stage4') && ...
            isfield(cached.stageDiagnostics.stage4,'radius')
        stageDiagnostics = cached.stageDiagnostics;
        return
    end
end

% Older caches pre-date the compact intermediate diagnostic payload. Build
% it once for the middle processed strake, then append it without changing
% any numerical protocol output.
fprintf(['Intermediate diagnostic data are absent from this older cache; ', ...
    'reconstructing one representative strake.\n'])
s1 = load(s1Cache,'regParams','bestFitParamsPreReg','iUse');
s2 = load(s2Cache,'JOINTBOUNDS');
s5 = load(s5Cache,'cloudStrakes','Z_BOTS','Z_TOPS','R0_BOTS','R0_TOPS', ...
    'THICKS','cacheConfig');
configuration = s5.cacheConfig;
[X,Y,Z] = loadRegisteredDiagnosticSample(fileparts(s5Cache), ...
    configuration.nperc,minZ,maxZ,s1.bestFitParamsPreReg, ...
    s1.regParams,s1.iUse);
bestFitParams = bestFitCone(X,Y,Z,minZ,tower,1e6);
[~,~,Z,T,R] = bestFitConeTransform(X,Y,Z,minZ,bestFitParams,false);
[~,boundaryRows] = selectCloudStrakes(s5.Z_BOTS,s5.Z_TOPS,minZ,maxZ);
if size(s2.JOINTBOUNDS,1) ~= numel(boundaryRows)
    error('generateProtocolDiagnostics:InvalidJointMapping', ...
        'Cached joint bounds do not match nominal boundaries in the selected height range.')
end
jointBounds = nan(numel(s5.Z_BOTS)+1,2);
jointBounds(boundaryRows,:) = s2.JOINTBOUNDS*1000;
strake = s5.cloudStrakes(ceil(numel(s5.cloudStrakes)/2));
bottomExtent = jointBounds(strake,2);
topExtent = jointBounds(strake+1,1);
selection = Z >= bottomExtent & Z <= topExtent;
thetaBefore = T(selection); zBefore = Z(selection); radiusBefore = R(selection);

expectedBeadWidth = configuration.upstream.expectedBeadWidth;
thetaBinSize = expectedBeadWidth*configuration.windowWidthSF*1000;
zBinSize = thetaBinSize;
if topExtent-bottomExtent < zBinSize
    zBinSize = topExtent-bottomExtent-5;
end
[thetaAfter,zAfter,radiusAfter] = S3_SystematicAndRandomOutlierRemoval( ...
    thetaBefore,zBefore,radiusBefore, ...
    (s5.R0_BOTS(strake)+s5.R0_TOPS(strake))/2,thetaBinSize,zBinSize, ...
    configuration.cutOff,false);

beforeIndices = displayIndices(numel(thetaBefore),50000);
afterIndices = displayIndices(numel(thetaAfter),50000);
stageDiagnostics.stage3.strake = strake;
stageDiagnostics.stage3.thetaBefore = thetaBefore(beforeIndices);
stageDiagnostics.stage3.zBefore = zBefore(beforeIndices);
stageDiagnostics.stage3.radiusBefore = radiusBefore(beforeIndices);
stageDiagnostics.stage3.thetaAfter = thetaAfter(afterIndices);
stageDiagnostics.stage3.zAfter = zAfter(afterIndices);
stageDiagnostics.stage3.radiusAfter = radiusAfter(afterIndices);
stageDiagnostics.stage3.inputCount = numel(thetaBefore);
stageDiagnostics.stage3.retainedCount = numel(thetaAfter);

[thetaGrid,zGrid,radiusGrid] = S4_GriddedOuterShellSurfaceReconstruction( ...
    thetaAfter,radiusAfter,zAfter,bottomExtent,topExtent,s5.R0_BOTS, ...
    s5.R0_TOPS,s5.Z_BOTS,s5.Z_TOPS,s5.THICKS,strake, ...
    configuration.searchRadiusGridSpacingSF,configuration.p, ...
    configuration.smoothingFilterStd,false);
rows = displayIndices(size(radiusGrid,1),300);
columns = displayIndices(size(radiusGrid,2),480);
stageDiagnostics.stage4.strake = strake;
stageDiagnostics.stage4.theta = thetaGrid(rows,columns);
stageDiagnostics.stage4.z = zGrid(rows,columns);
stageDiagnostics.stage4.radius = radiusGrid(rows,columns);
save(s5Cache,'stageDiagnostics','-append')
end

function [X,Y,Z] = loadRegisteredDiagnosticSample(inputDir,fraction,minZ,maxZ, ...
    bestFitParameters,registrationParameters,iterations)
files = dir(fullfile(inputDir,'*.bin'));
Xparts = cell(numel(files),1); Yparts = Xparts; Zparts = Xparts;
for cloudIndex = 1:numel(files)
    fprintf('  Sampling intermediate cloud %g of %g.\n',cloudIndex,numel(files))
    [x,y,z] = readBinaryPointSample( ...
        fullfile(files(cloudIndex).folder,files(cloudIndex).name),fraction,300000);
    [x,y,z] = bestFitConeTransform(x,y,z,minZ,bestFitParameters,false);
    if iterations > 0
        [x,y,z] = applyRegistrationTransform( ...
            x,y,z,registrationParameters,cloudIndex,iterations);
    end
    selection = z >= minZ & z <= maxZ;
    Xparts{cloudIndex} = x(selection);
    Yparts{cloudIndex} = y(selection);
    Zparts{cloudIndex} = z(selection);
end
X = vertcat(Xparts{:}); Y = vertcat(Yparts{:}); Z = vertcat(Zparts{:});
end

function plotOutlierRemovalDiagnostic(inputDir,tower,data)
fprintf('Generating Stage 3 outlier-removal diagnostic.\n')
referenceRadius = median(data.radiusAfter,'omitnan');
beforeDeviation = data.radiusBefore-referenceRadius;
afterDeviation = data.radiusAfter-referenceRadius;
beforeFinite = beforeDeviation(isfinite(beforeDeviation));
afterFinite = afterDeviation(isfinite(afterDeviation));
finiteDeviation = abs([beforeFinite(:);afterFinite(:)]);
colourLimit = diagnosticColourLimit(finiteDeviation);

figureHandle = figure('Visible','off','Color','w','Position',[100 100 1400 600]);
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact','Padding','compact');
for panel = 1:2
    axesHandle = nexttile(layout);
    if panel == 1
        theta = data.thetaBefore; z = data.zBefore; deviation = beforeDeviation;
        panelTitle = sprintf('Input display sample (%s points)',formatCount(numel(theta)));
    else
        theta = data.thetaAfter; z = data.zAfter; deviation = afterDeviation;
        panelTitle = sprintf('Retained display sample (%s points)',formatCount(numel(theta)));
    end
    scatter(axesHandle,theta,z/1000,6,deviation,'filled','MarkerFaceAlpha',0.65)
    xlim(axesHandle,[0 2*pi]); clim(axesHandle,[-colourLimit colourLimit]);
    grid(axesHandle,'on'); box(axesHandle,'on'); colormap(axesHandle,turbo)
    xlabel(axesHandle,'$\theta$ [rad]','Interpreter','latex')
    ylabel(axesHandle,'Elevation $z$ [m]','Interpreter','latex')
    title(axesHandle,panelTitle,'Interpreter','latex')
    colourBar = colorbar(axesHandle);
    colourBar.Label.String = 'Radial deviation [mm]';
    styleColourBar(colourBar)
    styleAxes(axesHandle)
end
inputCount = diagnosticCount(data.inputCount,data.strake);
retainedCount = diagnosticCount(data.retainedCount,data.strake);
retained = 100*retainedCount/max(inputCount,1);
title(layout,sprintf('Stage 3 outlier removal: representative strake %g (%.1f\\%% retained)', ...
    data.strake,retained),'Interpreter','latex','Color','k')
diagnosticCaption(layout,{'Before and after outlier removal for one strake, using the same colour scale.', ...
    'Check for loss of coherent measured features as well as removal of scattered points.'})
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(3,tower))
end

function plotGriddedSurfaceDiagnostic(inputDir,tower,data)
fprintf('Generating Stage 4 gridded-surface diagnostic.\n')
deviation = data.radius-median(data.radius,2,'omitnan');
finiteDeviation = abs(deviation(isfinite(deviation)));
colourLimit = diagnosticColourLimit(finiteDeviation);

figureHandle = figure('Visible','off','Color','w','Position',[100 100 1400 600]);
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact','Padding','compact');
axesHandle = nexttile(layout);
surf(axesHandle,data.theta,data.z/1000,deviation,'EdgeColor','none')
view(axesHandle,2); axis(axesHandle,'tight'); xlim(axesHandle,[0 2*pi]);
clim(axesHandle,[-colourLimit colourLimit]); colormap(axesHandle,turbo)
grid(axesHandle,'on'); box(axesHandle,'on')
xlabel(axesHandle,'$\theta$ [rad]','Interpreter','latex')
ylabel(axesHandle,'Local $z$ [m]','Interpreter','latex')
title(axesHandle,'Reconstructed radial deviation','Interpreter','latex')
colourBar = colorbar(axesHandle);
colourBar.Label.String = 'Radial deviation from row median [mm]';
styleColourBar(colourBar); styleAxes(axesHandle)

axesHandle = nexttile(layout); hold(axesHandle,'on')
profileRows = unique(round(linspace(1,size(deviation,1),3)));
profileColours = lines(numel(profileRows));
profileLabels = cell(numel(profileRows),1);
for profileIndex = 1:numel(profileRows)
    row = profileRows(profileIndex);
    plot(axesHandle,data.theta(row,:),deviation(row,:),'LineWidth',1.25, ...
        'Color',profileColours(profileIndex,:))
    profileLabels{profileIndex} = sprintf('$z = %.3f$ m',data.z(row,1)/1000);
end
xlim(axesHandle,[0 2*pi]); ylim(axesHandle,[-colourLimit colourLimit]);
grid(axesHandle,'on'); box(axesHandle,'on')
xlabel(axesHandle,'$\theta$ [rad]','Interpreter','latex')
ylabel(axesHandle,'Radial deviation from row median [mm]','Interpreter','latex')
title(axesHandle,'Representative circumferential profiles','Interpreter','latex')
legendHandle = legend(axesHandle,profileLabels,'Location','best', ...
    'Interpreter','latex');
set(legendHandle,'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25])
styleAxes(axesHandle)
title(layout,sprintf('Stage 4 gridded outer surface: representative strake %g', ...
    data.strake),'Interpreter','latex','Color','k')
diagnosticCaption(layout,{'Only nodes with neighbours inside the fixed search radius are reconstructed.', ...
    'Blank areas have no local support and remain NaN after smoothing.'})
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(4,tower))
end

function plotSurfaceDiagnostic(inputDir,tower,s5Cache)
fprintf('Generating Stage 5 reconstructed-surface diagnostic.\n')
s5 = load(s5Cache,'cloudStrakes','CANS');
allDeviation = [];
for strake = s5.cloudStrakes(:)'
    radii = s5.CANS.(sprintf('S%g',strake)).RBIG;
    deviation = radii-median(radii,2,'omitnan');
    allDeviation = [allDeviation; deviation(:)]; %#ok<AGROW>
end
colourLimit = diagnosticColourLimit(allDeviation);

figureHandle = figure('Visible','off','Color','w','Position',[100 100 1250 720]);
axesHandle = axes(figureHandle); hold(axesHandle,'on')
for strake = s5.cloudStrakes(:)'
    can = s5.CANS.(sprintf('S%g',strake));
    deviation = can.RBIG-median(can.RBIG,2,'omitnan');
    rows = displayIndices(size(deviation,1),300);
    columns = displayIndices(size(deviation,2),480);
    surf(axesHandle,can.THETg(rows,columns),can.Z3D(rows,columns)/1000, ...
        deviation(rows,columns),'EdgeColor','none');
end
view(axesHandle,2); axis(axesHandle,'tight'); box(axesHandle,'on')
xlim(axesHandle,[0 2*pi]); clim(axesHandle,[-colourLimit colourLimit]);
colormap(axesHandle,turbo)
xlabel(axesHandle,'$\theta$ [rad]','Interpreter','latex');
ylabel(axesHandle,'Elevation $z$ [m]','Interpreter','latex')
title(axesHandle,'Stage 5 reconstructed surface: circumferential radial deviation', ...
    'Interpreter','latex')
diagnosticCaption(axesHandle,{'Colour shows radial deviation from each row median, not measurement density.', ...
    'Blank areas remain unsupported; compare with the measured point-cloud figure.'})
colourBar = colorbar(axesHandle);
colourBar.Label.String = 'Radial deviation from row median [mm]';
colourBar.Label.Interpreter = 'latex';
colourBar.TickLabelInterpreter = 'latex';
colourBar.FontSize = 12;
colourBar.Color = 'k';
styleAxes(axesHandle)
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(5,tower))
end

function plotPointCloudDiagnostic(inputDir,minZ,maxZ,s1Cache,s5Cache,stem)
fprintf('Generating Stage 5 measured point-cloud diagnostic.\n')
s1 = load(s1Cache,'regParams','bestFitParamsPreReg','iUse');
s5 = load(s5Cache,'cacheConfig');
% Use the registration selected for this surface, including any override.
if isfield(s5.cacheConfig.upstream,'registrationIteration')
    s1.iUse = s5.cacheConfig.upstream.registrationIteration;
end
files = dir(fullfile(inputDir,'*.bin'));
if isempty(files)
    error('generateProtocolDiagnostics:NoClouds','No binary point clouds in %s.',inputDir)
end
colours = lines(numel(files));
figureHandle = figure('Visible','off','WindowStyle','normal','Color','w', ...
    'Position',[100 100 1400 850]);
layout = tiledlayout(figureHandle,1,2,'TileSpacing','compact','Padding','compact');
axes3D = nexttile(layout); hold(axes3D,'on')
axesUnwrapped = nexttile(layout); hold(axesUnwrapped,'on')
handles = gobjects(numel(files),1);
labels = cell(numel(files),1);
for scan = 1:numel(files)
    % Deterministic display sample of measured XYZ, before outlier removal
    % and reconstruction. Keep each scan separate to preserve its identity.
    [x,y,z] = readBinaryPointSample(fullfile(files(scan).folder,files(scan).name),1,50000);
    [x,y,z] = bestFitConeTransform(x,y,z,minZ,s1.bestFitParamsPreReg,false);
    if s1.iUse > 0
        [x,y,z] = applyRegistrationTransform(x,y,z,s1.regParams,scan,s1.iUse);
    end
    selection = z >= minZ & z <= maxZ;
    x = x(selection); y = y(selection); z = z(selection);
    [~,labels{scan}] = fileparts(files(scan).name);
    handles(scan) = scatter3(axes3D,x,y,z,3,colours(scan,:),'filled', ...
        'DisplayName',labels{scan});
    scatter(axesUnwrapped,mod(atan2(y,x),2*pi),z,3,colours(scan,:),'filled')
    fprintf('  %s: %g displayed measured points.\n',labels{scan},numel(z))
end
axis(axes3D,'equal'); axis(axes3D,'tight'); view(axes3D,38,24)
xlabel(axes3D,'$x$ [m]','Interpreter','latex')
ylabel(axes3D,'$y$ [m]','Interpreter','latex')
zlabel(axes3D,'Elevation $z$ [m]','Interpreter','latex')
title(axes3D,'Registered measured points','Interpreter','latex')
xlim(axesUnwrapped,[0 2*pi]); ylim(axesUnwrapped,[minZ maxZ])
xlabel(axesUnwrapped,'$\theta$ [rad]','Interpreter','latex')
ylabel(axesUnwrapped,'Elevation $z$ [m]','Interpreter','latex')
title(axesUnwrapped,'Measured circumferential coverage','Interpreter','latex')
for axesHandle = [axes3D axesUnwrapped]
    grid(axesHandle,'on'); box(axesHandle,'on'); styleAxes(axesHandle)
end
legendHandle = legend(axes3D,handles,labels,'Interpreter','none', ...
    'Orientation','horizontal','NumColumns',min(4,numel(files)));
legendHandle.Layout.Tile = 'south';
set(legendHandle,'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25])
title(layout,'Stage 5: measured point cloud by scan position', ...
    'Interpreter','latex','Color','k')
diagnosticCaption(layout,{'Registered measurements before outlier removal, sampled to at most 50,000 points per scan.', ...
    'Colours identify scan files. Use the unwrapped view to inspect gaps hidden on the far side in 3D.'})
saveDiagnosticFigure(figureHandle,inputDir,stem)
end

function plotMeshDiagnostic(inputDir,tower,s6Cache)
fprintf('Generating Stage 6 mesh diagnostic.\n')
s6 = load(s6Cache,'X_MESH','Y_MESH','Z_MESH','R_MESH','cloudStrakes','N2_MESH');
deviation = s6.R_MESH-median(s6.R_MESH,2,'omitnan');
colourLimit = diagnosticColourLimit(deviation);
columns = displayIndices(size(deviation,2),720);
figureHandle = figure('Visible','off','Color','w','Position',[100 100 720 1000]);
axesHandle = axes(figureHandle); hold(axesHandle,'on')
firstRow = 1;
for strake = s6.cloudStrakes(:)'
    rows = firstRow-1+displayIndices(s6.N2_MESH(strake),500);
    firstRow = firstRow+s6.N2_MESH(strake);
    % Close the circumferential seam for display; keep strakes separate.
    xPlot = [s6.X_MESH(rows,columns),s6.X_MESH(rows,1)];
    yPlot = [s6.Y_MESH(rows,columns),s6.Y_MESH(rows,1)];
    zPlot = [s6.Z_MESH(rows,columns),s6.Z_MESH(rows,1)];
    deviationPlot = [deviation(rows,columns),deviation(rows,1)];

    surf(axesHandle,xPlot/1000,yPlot/1000,zPlot/1000,deviationPlot,'EdgeColor','none');
end
axis(axesHandle,'equal'); axis(axesHandle,'tight'); box(axesHandle,'on');
grid(axesHandle,'on'); view(axesHandle,38,24)
clim(axesHandle,[-colourLimit colourLimit]); colormap(axesHandle,turbo)
xlabel(axesHandle,'$x$ [m]','Interpreter','latex');
ylabel(axesHandle,'$y$ [m]','Interpreter','latex');
zlabel(axesHandle,'$z$ [m]','Interpreter','latex')
title(axesHandle,'Stage 6 projected shell mesh','Interpreter','latex')
diagnosticCaption(axesHandle,{'The reconstructed surface is interpolated onto the target mesh.', ...
    'Missing regions remain NaN; a mesh with holes is not a complete shell model.'})
colourBar = colorbar(axesHandle);
colourBar.Label.String = 'Radial deviation from row median [mm]';
colourBar.Label.Interpreter = 'latex';
colourBar.TickLabelInterpreter = 'latex';
colourBar.FontSize = 12;
colourBar.Color = 'k';
styleAxes(axesHandle)
saveDiagnosticFigure(figureHandle,inputDir,diagnosticStem(6,tower))
end

function diagnosticCaption(target,lines)
subtitle(target,lines,'Interpreter','none','FontSize',11,'Color',[0.25 0.25 0.25])
end

function limit = diagnosticColourLimit(values)
values = abs(values(isfinite(values)));
if isempty(values)
    limit = 1; % Allow an entirely unsupported surface to display as empty.
else
    limit = max(prctile(values,99),eps);
end
end

function indices = displayIndices(count,maximumCount)
indices = unique(round(linspace(1,count,min(count,maximumCount))));
end

function count = diagnosticCount(values,strake)
if isscalar(values)
    count = values;
else
    count = values(strake);
end
end

function value = formatCount(count)
value = sprintf('%g',count);
end

function styleColourBar(colourBar)
colourBar.Label.Interpreter = 'latex';
colourBar.TickLabelInterpreter = 'latex';
colourBar.FontSize = 12;
colourBar.Color = 'k';
end

function saveDiagnosticFigure(figureHandle,inputDir,stem)
set(figureHandle,'Color','w','InvertHardcopy','off')
% Save the interactive figure before raster rendering so the MATLAB output
% is preserved even if a machine's graphics service cannot export a PNG.
% FIG files must open visibly when selected in MATLAB or Explorer.
figureHandle.Visible = 'on';
savefig(figureHandle,fullfile(inputDir,[stem,'.fig']))
exportgraphics(figureHandle,fullfile(inputDir,[stem,'.png']), ...
    'Resolution',300,'BackgroundColor','white')
close(figureHandle)
end

function styleAxes(axesHandle)
set(axesHandle,'Color','w','XColor','k','YColor','k','ZColor','k', ...
    'GridColor',[0.75 0.75 0.75],'MinorGridColor',[0.85 0.85 0.85], ...
    'TickLabelInterpreter','latex','FontSize',12,'LineWidth',0.8)
axesHandle.Title.Color = 'k';
axesHandle.XLabel.Color = 'k';
axesHandle.YLabel.Color = 'k';
axesHandle.ZLabel.Color = 'k';
end
