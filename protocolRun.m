% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Run the code to go through each stage of the protocol. Once a stage has
% been completed, the relevant data is saved to a .mat file in the same
% directory as the point clouds. The protocol therfore does not need to be
% run from the start each time for a dataset.

clear;
clc;
close all;

%% USER SETTINGS
% Set these for the dataset before running. The joint count comes from the
% nominal geometry, as do flange types and the eligible strake selection.

% All stages: dataset, elevation range and restart
tower = 'F3'; % Input_<tower> folder and S0_Data_<tower>.m geometry file
minZ = -Inf; % [m] -Inf uses the lowest measured elevation
maxZ = Inf; % [m] Inf uses the highest measured elevation
startFromStage = 1; % 1--7; stages 3--5 run together. 7 loads completed outputs.
randomSeed = 1729; % Random sampling and optimisation, all stages

% Diagnostics: all stages
plotting = false; % Detailed figures during processing; can be slow
generateDiagnostics = true; % Standard FIG/PNG outputs for stages 1--6
forceDiagnostics = false; % Rebuild existing standard figures
units = 'm'; % Stage 1 plot labels; input coordinates must be in metres

% Stage 1: enhanced registration
registrationFraction = 0.01; % Fraction of each scan retained
numIter = 5; % Registration iterations
cutoffReg = 0.02; % [m] Maximum distance for registration correspondences
registrationIteration = []; % [] selects the best iteration; 1--numIter overrides it for stages 2--6

% Stage 2: radial profile and joint detection
profileFraction = 0.1; % Fraction of each scan retained
expectedBeadWidth = 0.08; % [m] Expected weld bead width; also sets Stage 3 window sizes
numCircWindows = 100; % Circumferential windows for radial correction
subIntervals = 10; % Vertical subintervals per profile window; at least 2
jointNominalTolerance = 0.10; % [m] Maximum measured/nominal joint offset; validation only
jointSearchRadius = 0.10; % [m] Peak search on either side of each nominal joint, clipped at adjacent midpoints
trimFailedEndJoints = true; % Stage 2: allow incomplete end profiles and missing end junctions to be excluded
retainLongestProfileRun = true; % Stage 2: with end trimming enabled, keep the longest continuous profile; otherwise interior gaps stop
retainJointRuns = true; % Stage 2: keep separate runs of detected joints; exclude strakes with a missing boundary
profileEndTrimDistance = 0.5; % [m] Patchy-end allowance when retainLongestProfileRun is false; 0 disables
% Joint count and flange types come from Z_BOTS and IS_FLANGE in S0_Data_<tower>.m.
flangePadding = 10; % [mm] Extend fitted flange-joint bounds by this amount
weldPadding = 5; % [mm] Extend fitted weld bounds by this amount

% Stages 3--6: surface selection and sampling
% cloudStrakes is automatic: segments with two internal junctions inside the height range.
% Tower-end segments and segments cut by the height limits are excluded.
surfaceFraction = 0.02; % Fraction of each scan used for stages 3--5

% Stage 3: outlier removal
cutOff = 1; % [mm] Point-to-plane outlier threshold
windowWidthSF = 3; % Window size / expected bead width

% Stage 4: surface reconstruction
p = 1; % Inverse-distance weighting exponent
searchRadiusGridSpacingSF = 4; % Fixed search radius / estimated point spacing
smoothingFilterStd = 10; % Gaussian standard deviation, in grid cells

% Stage 5 uses the Stage 2 joint bounds and the Stage 4 surface.
% Stage 6 derives mesh spacing from nominal radii, heights and thicknesses.
% Neither has additional user settings in this script.

% Cache versions: change after numerical implementation changes
implementationVersion = 2; % All stages
reconstructionVersion = 6; % Stages 3--6: explicit result for an unusable cone fit

%% SETUP
rng(randomSeed,'twister')
repoRoot = fileparts(mfilename('fullpath'));
inputDir = fullfile(repoRoot,['Input_',tower]);
if ~isempty(registrationIteration)
    validateattributes(registrationIteration,{'numeric'}, ...
        {'scalar','integer','>=',1,'<=',numIter},mfilename,'registrationIteration')
end

%% Adding folders to path
addpath(genpath(fullfile(repoRoot,'Common')))
addpath(genpath(inputDir))
addpath(genpath(fullfile(repoRoot,'Stage_1')))
addpath(genpath(fullfile(repoRoot,'Stage_2')))
addpath(genpath(fullfile(repoRoot,'Stage_3')))
addpath(genpath(fullfile(repoRoot,'Stage_4')))
addpath(genpath(fullfile(repoRoot,'Stage_5')))
addpath(genpath(fullfile(repoRoot,'Stage_6')))

%% Nominal input data
clear IS_FLANGE
run(fullfile(inputDir,['S0_Data_',tower,'.m']))
if ~exist('IS_FLANGE','var')
    error('protocolRun:MissingFlangeTypes', ...
        'Add IS_FLANGE to S0_Data_%s.m: one true/false value per nominal segment. Segment names and thicknesses are not used to infer flange types.',tower)
end
if ~((islogical(IS_FLANGE) || isnumeric(IS_FLANGE)) && isreal(IS_FLANGE) && ...
        isvector(IS_FLANGE) && numel(IS_FLANGE) == numel(Z_BOTS) && ...
        all(IS_FLANGE(:) == 0 | IS_FLANGE(:) == 1))
    error('protocolRun:InvalidFlangeTypes', ...
        'IS_FLANGE in S0_Data_%s.m must be a logical or 0/1 vector with one entry per Z_BOTS segment.',tower)
end
isFlange = logical(IS_FLANGE(:)');

%% Create .bin files
convertPTS2BIN(tower) % Creates .bin files if no .bin files exist in dataset directory
inputSignature = pointCloudFileSignature(inputDir);
[minZ,maxZ] = resolveHeightRange(inputDir,minZ,maxZ); % Resolve before fitting or cache checks.
cloudStrakes = selectCloudStrakes(Z_BOTS,Z_TOPS,minZ,maxZ);
fprintf('Selected nominal strakes: %s\n',mat2str(cloudStrakes))
fprintf('Tower-end segments and segments without two in-range junctions are excluded.\n')

%% STAGE 1: ENHANCED REGISTRATION
nperc = registrationFraction;
s1Cache = fullfile(inputDir,['S1_Registration_',tower,'.mat']);
stage1Config = struct('stage',1,'implementationVersion',implementationVersion, ...
    'inputSignature',inputSignature,'randomSeed',randomSeed,'minZ',minZ, ...
    'maxZ',maxZ,'nperc',nperc,'numIter',numIter,'cutoffReg',cutoffReg);

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
if startFromStage > 1 && isCompatibleCache(s1Cache,stage1Config)
    disp('STAGE 1: ENHANCED REGISTRATION - Loading')
    loaded = load(s1Cache,'nClouds','regParams','bestFitParamsPreReg','iUse');
    nClouds = loaded.nClouds; regParams = loaded.regParams;
    bestFitParamsPreReg = loaded.bestFitParamsPreReg; iUse = loaded.iUse;
else
    disp('STAGE 1: ENHANCED REGISTRATION - Start')
    [regParams, bestFitParamsPreReg, iUse, nClouds] = S1_EnhancedRegistration(minZ, maxZ, units, nperc, tower, numIter, cutoffReg, plotting); % Doing registration
    cacheConfig = stage1Config;
    save(s1Cache,'nClouds','regParams','bestFitParamsPreReg','iUse','cacheConfig','-v7.3');
    disp('STAGE 1: ENHANCED REGISTRATION - Complete')
end

%% STAGE 2: CAN SEGMENTATION
% Internal nominal junctions within the selected height range [m].
nominalJointZ = Z_BOTS(2:end)/1000;
inRange = nominalJointZ > minZ & nominalJointZ < maxZ;
% Only the shared boundary of two flange segments is a flange connection.
nominalFlangeJoint = isFlange(1:end-1) & isFlange(2:end);
nominalFlangeJoint = nominalFlangeJoint(inRange);
nominalJointZ = nominalJointZ(inRange);
nJoints = numel(nominalJointZ);
fprintf('Stage 2: %g nominal internal junctions in the selected height range.\n',nJoints)
nperc = profileFraction;
if ~isempty(registrationIteration)
    iUse = registrationIteration;
end
s2Cache = fullfile(inputDir,['S2_CanSegmentation_',tower,'.mat']);
validateattributes(profileEndTrimDistance,{'numeric'},{'scalar','real','finite','nonnegative'})
profileInteriorRange = [-Inf Inf];
if profileEndTrimDistance > 0
    profileInteriorRange = [Z_BOTS(1)/1000+profileEndTrimDistance, Z_TOPS(end)/1000-profileEndTrimDistance];
end
stage2Config = struct('stage',2,'implementationVersion',implementationVersion, ...
    'profileVersion',7,'profileInteriorRange',profileInteriorRange, ...
    'retainJointRuns',retainJointRuns, ...
    'retainLongestProfileRun',retainLongestProfileRun, ...
    'upstream',stage1Config,'nperc',nperc,'expectedBeadWidth',expectedBeadWidth, ...
    'numCircWindows',numCircWindows,'subIntervals',subIntervals,'nJoints',nJoints, ...
    'nominalJointZ',nominalJointZ,'nominalFlangeJoint',nominalFlangeJoint, ...
    'jointSearchRadius',jointSearchRadius,'trimFailedEndJoints',trimFailedEndJoints, ...
    'flangePadding',flangePadding,'weldPadding',weldPadding);
if ~isempty(registrationIteration)
    stage2Config.registrationIteration = registrationIteration;
end

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
if startFromStage > 2 && isCompatibleCache(s2Cache,stage2Config)
    disp('STAGE 2: CAN SEGMENTATION - Loading')
    loaded = load(s2Cache,'ZMEDIANS','RMEDIANS','jointCoordinates','JOINTWIDTH','JOINTBOUNDS');
    ZMEDIANS = loaded.ZMEDIANS; RMEDIANS = loaded.RMEDIANS;
    jointCoordinates = loaded.jointCoordinates; JOINTWIDTH = loaded.JOINTWIDTH;
    JOINTBOUNDS = loaded.JOINTBOUNDS;
else

    disp('STAGE 2: CAN SEGMENTATION - Start')

    % Global radial profile generation 
    disp('STAGE 2: Global radial profile generation ')
    [ZMEDIANS, RMEDIANS] = S2_GlobalRadialProfileGeneration(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower, expectedBeadWidth, numCircWindows, subIntervals, plotting, trimFailedEndJoints, profileInteriorRange, retainLongestProfileRun);
    
    % Joint ROI identification
    disp('STAGE 2: Joint ROI identification')
    [jointCoordinates,jointMatches] = S2_JointROIIdentification( ...
        ZMEDIANS,RMEDIANS,nJoints,expectedBeadWidth,plotting,nominalJointZ,jointSearchRadius);
    if trimFailedEndJoints || retainJointRuns
        outsideProfile = nominalJointZ < ZMEDIANS(1) | nominalJointZ > ZMEDIANS(end);
        jointCoordinates(outsideProfile) = NaN;
        jointMatches{outsideProfile,{'Detected_m','Offset_m','PeakStrength_m'}} = NaN;
    end
    % Preserve the measured result even if validation rejects it.
    retainedJoints = retainedJointRun(jointCoordinates,trimFailedEndJoints,retainJointRuns,nominalFlangeJoint,nominalJointZ);
    detectionConfig = stage2Config;
    save(fullfile(inputDir,['S2_JointDetection_',tower,'.mat']), ...
        'ZMEDIANS','RMEDIANS','jointCoordinates','nominalJointZ','jointMatches','retainedJoints','detectionConfig','-v7.3')
    if generateDiagnostics
        try
            generateProtocolDiagnostics(tower,minZ,maxZ,true,2)
        catch diagnosticError
            warning('protocolRun:Stage2DiagnosticGenerationFailed', ...
                'Stage 2 detection data are saved, but plotting failed: %s',diagnosticError.message)
        end
    end
    if nnz(retainedJoints) < 2
        error('protocolRun:NoCompleteStrakes','At least two valid junctions are needed to bound a strake.')
    end
    validateJointDetection(jointCoordinates(retainedJoints),nominalJointZ(retainedJoints),jointNominalTolerance)
    flangeTransitionJointIndex = matchFlangeJoints( ...
        jointCoordinates(retainedJoints),nominalJointZ(retainedJoints),nominalFlangeJoint(retainedJoints));

    % Joint segmentation
    disp('STAGE 2: Joint segmentation')
    [widths,bounds] = S2_JointSegmentation(ZMEDIANS, RMEDIANS, jointCoordinates(retainedJoints), flangeTransitionJointIndex, flangePadding, weldPadding, plotting);
    % Keep every nominal row, including excluded ends, in saved results.
    JOINTWIDTH = nan(nJoints,1); JOINTBOUNDS = nan(nJoints,2);
    JOINTWIDTH(retainedJoints) = widths; JOINTBOUNDS(retainedJoints,:) = bounds;
    retainedIndices = find(retainedJoints);
    flangeTransitionJointIndex = retainedIndices(flangeTransitionJointIndex);
    cacheConfig = stage2Config;
    save(s2Cache,'ZMEDIANS','RMEDIANS','jointCoordinates','JOINTWIDTH','JOINTBOUNDS', ...
        'flangeTransitionJointIndex','nominalJointZ','nominalFlangeJoint','jointMatches','retainedJoints','cacheConfig','-v7.3');

    disp('STAGE 2: CAN SEGMENTATION - Complete')
end

%% Preparing data for looping through cans
% Stage 2 plots must not depend on later stages completing successfully.
if generateDiagnostics
    try
        generateProtocolDiagnostics(tower,minZ,maxZ,forceDiagnostics,2)
    catch diagnosticError
        warning('protocolRun:Stage2DiagnosticGenerationFailed', ...
            'Stage 2 data are saved, but plotting failed: %s',diagnosticError.message)
    end
end
% Check cached results too, before loading the point clouds or surfaces.
retainedJoints = retainedJointRun(jointCoordinates,trimFailedEndJoints,retainJointRuns,nominalFlangeJoint,nominalJointZ);
validateJointDetection(jointCoordinates(retainedJoints),nominalJointZ(retainedJoints),jointNominalTolerance)
[~,boundaryRows] = selectCloudStrakes(Z_BOTS,Z_TOPS,minZ,maxZ);
availableBoundaries = false(numel(Z_BOTS)+1,1);
availableBoundaries(boundaryRows(retainedJoints)) = true;
cloudStrakes = find(availableBoundaries(1:end-1) & availableBoundaries(2:end))';
if isempty(cloudStrakes)
    error('protocolRun:NoCompleteStrakes','No strake has two retained junctions.')
end
if any(~retainedJoints)
    fprintf('Stage 2: excluded junctions at nominal elevations %s m (missing peak, flange neighbours or isolated joint).\n',mat2str(nominalJointZ(~retainedJoints)))
    fprintf('Retained nominal strakes: %s; detected bounding centres %.3f--%.3f m.\n', ...
        mat2str(cloudStrakes),jointCoordinates(find(retainedJoints,1)),jointCoordinates(find(retainedJoints,1,'last')))
end
runStarts = find(retainedJoints & [true ~retainedJoints(1:end-1)]);
runEnds = find(retainedJoints & [~retainedJoints(2:end) true]);
for runIndex = 1:numel(runStarts)
    fprintf('Stage 2: retained joint run %.3f--%.3f m; only strakes with both boundaries are reconstructed.\n', ...
        jointCoordinates(runStarts(runIndex)),jointCoordinates(runEnds(runIndex)))
end
strakeJointBounds = validateStrakeSelection( ...
    cloudStrakes,Z_BOTS,Z_TOPS,JOINTBOUNDS,minZ,maxZ,jointNominalTolerance,retainedJoints,jointCoordinates);
nperc = surfaceFraction;
s5Cache = fullfile(inputDir,['S5_Surface_',tower,'.mat']);
stage5Config = struct('stage',[3,4,5],'implementationVersion',implementationVersion, ...
    'reconstructionVersion',reconstructionVersion, ...
    'upstream',stage2Config,'randomSeed',randomSeed,'cloudStrakes',cloudStrakes, ...
    'nperc',nperc,'cutOff',cutOff,'windowWidthSF',windowWidthSF,'p',p, ...
    'searchRadiusGridSpacingSF',searchRadiusGridSpacingSF, ...
    'smoothingFilterStd',smoothingFilterStd);
loadStage5 = startFromStage > 5 && isCompatibleCache(s5Cache,stage5Config);

% If a later start stage is requested and the relevant .mat file exists, no need to load point cloud
if loadStage5
    disp('')
else
    
    % Loading data
    [X, Y, Z, ~, ~] = postRegDataLoader(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower);
    
    % Performing best fit cone
    [bestFitParams] = bestFitCone(X, Y, Z, minZ, tower, 1e6);
    [X, Y, Z, T, R] = bestFitConeTransform(X, Y, Z, minZ, bestFitParams);
    
    % Loading tower details
    run(fullfile(inputDir,['S0_Data_',tower]))
    
    
    % Applying unit conversion of point cloud to mm
    X = X*1000; Y = Y*1000; Z = Z*1000; R = R*1000;
    JOINTBOUNDS = strakeJointBounds*1e3;
    N1 = zeros(size(Z_BOTS));
    N2 = N1;
    diagnosticStrake = cloudStrakes(ceil(numel(cloudStrakes)/2));
    stageDiagnostics = struct();
    stageDiagnostics.stage3.strake = diagnosticStrake;
    stageDiagnostics.stage3.inputCount = nan(size(Z_BOTS));
    stageDiagnostics.stage3.retainedCount = nan(size(Z_BOTS));
    stageDiagnostics.stage4.strake = diagnosticStrake;
end

%% Looping through cans and performing stages 3, 4 and 5 for each can

% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stages
if loadStage5
    disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Loading')
    disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Loading')
    disp('STAGE 5: SURFACE INPAINTING - Loading')
    loaded = load(s5Cache,'cloudStrakes','CANS','N1','N2','Z_BOTS','Z_TOPS','R0_BOTS','R0_TOPS','THICKS');
    cloudStrakes = loaded.cloudStrakes; CANS = loaded.CANS;
    N1 = loaded.N1; N2 = loaded.N2; Z_BOTS = loaded.Z_BOTS;
    Z_TOPS = loaded.Z_TOPS; R0_BOTS = loaded.R0_BOTS;
    R0_TOPS = loaded.R0_TOPS; THICKS = loaded.THICKS;
else
    for S = 1:length(Z_BOTS)
        disp(['Creating Strake ', num2str(S)])
    
        if ismember(S, cloudStrakes)
            TT = T;
            RR = R;
            ZZ = Z;
            XX = X;
            YY = Y;
           
            botExtent = JOINTBOUNDS(S,2); topExtent = JOINTBOUNDS(S+1,1); % Top and bottom of can
            JOINTCUTOFFTOP = ((JOINTBOUNDS(S+1,1)+JOINTBOUNDS(S+1,2))/2-JOINTBOUNDS(S+1,1)); % Distance from top of can to middle of above joint
            JOINTCUTOFFBOT = (JOINTBOUNDS(S,2)-((JOINTBOUNDS(S,1)+JOINTBOUNDS(S,2))/2)); % Distance from bottom of can to middle of below joint
    
            % Identifying can region
            sel = Z < botExtent | Z > topExtent;
            TT(sel) = [];
            ZZ(sel) = [];
            RR(sel) = [];
            XX(sel) = [];
            YY(sel) = [];
    
            % STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL
            disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Start')
            tbinSize = expectedBeadWidth*windowWidthSF*1e3; zbinSize = tbinSize; % Window sizes for outlier removal
            if topExtent-botExtent < zbinSize; zbinSize = topExtent-botExtent-5; end % Accounting for if can is smaller than zbinSize
            [TT_REM, ZZ_REM, RR_REM, ~, ~, ~, ~, ~, ~, ~] = S3_SystematicAndRandomOutlierRemoval(TT, ZZ, RR, (R0_BOTS(S) + R0_TOPS(S))/2, tbinSize, zbinSize, cutOff, plotting);
            stageDiagnostics.stage3.inputCount(S) = numel(TT);
            stageDiagnostics.stage3.retainedCount(S) = numel(TT_REM);
            if S == diagnosticStrake
                beforeIndices = unique(round(linspace(1,numel(TT),min(numel(TT),50000))));
                afterIndices = unique(round(linspace(1,numel(TT_REM),min(numel(TT_REM),50000))));
                stageDiagnostics.stage3.thetaBefore = TT(beforeIndices);
                stageDiagnostics.stage3.zBefore = ZZ(beforeIndices);
                stageDiagnostics.stage3.radiusBefore = RR(beforeIndices);
                stageDiagnostics.stage3.thetaAfter = TT_REM(afterIndices);
                stageDiagnostics.stage3.zAfter = ZZ_REM(afterIndices);
                stageDiagnostics.stage3.radiusAfter = RR_REM(afterIndices);
            end
            if plotting
                title(sprintf('Outlier removal of strake at %g - %g m',Z_BOTS(S)/1000, Z_TOPS(S)/1000),'interpreter','latex')
            end
            disp('STAGE 3: SYSTEMATIC AND RANDOM OUTLIER REMOVAL - Complete')

            % STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION
            disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Start')
            [THETg, Zg, RBIG, reconstructionFailure] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd, plotting);
            if ~isempty(reconstructionFailure)
                warning('protocolRun:SkippedStrake', ...
                    'Stage 4: skipping strake %g (%.3f--%.3f m): %s', ...
                    S,botExtent/1000,topExtent/1000,reconstructionFailure.reason)
                stageDiagnostics.stage4.skippedStrakes(S) = reconstructionFailure;
                cloudStrakes(cloudStrakes == S) = [];
                N1(S) = 0; N2(S) = 0;
                CANS.(sprintf('S%g',S)) = struct('THETg',[],'Zg',[],'RBIG',[], ...
                    'X3D',[],'Y3D',[],'Z3D',[],'N1',0,'N2',0);
                continue
            end
            if S == diagnosticStrake
                diagnosticRows = unique(round(linspace(1,size(RBIG,1),min(size(RBIG,1),300))));
                diagnosticColumns = unique(round(linspace(1,size(RBIG,2),min(size(RBIG,2),480))));
                stageDiagnostics.stage4.theta = THETg(diagnosticRows,diagnosticColumns);
                stageDiagnostics.stage4.z = Zg(diagnosticRows,diagnosticColumns);
                stageDiagnostics.stage4.radius = RBIG(diagnosticRows,diagnosticColumns);
            end
            disp('STAGE 4: GRIDDED OUTER SURFACE RECONSTRUCTION - Complete')

            % STAGE 5: SURFACE INPAINTING
            disp('STAGE 5: SURFACE INPAINTING - Start')
            [THETg, Zg, RBIG] = S5_SurfaceInpainting(THETg, Zg, RBIG, JOINTCUTOFFBOT, JOINTCUTOFFTOP);
            disp('STAGE 5: SURFACE INPAINTING - Complete')
    
            if plotting
                figure;
                surf(THETg,Zg,RBIG)
                grid on
                xlabel('$\theta$ [rad]','interpreter','latex')
                ylabel('$z$ [mm]','interpreter','latex')
                zlabel('$\rho$ [mm]','interpreter','latex')
                title(sprintf('Surface inpainting of strake at %g - %g m',Z_BOTS(S)/1000,Z_TOPS(S)/1000),'interpreter','latex')
            end
    
            % Computing cartesian coordinates of can
            X3D = RBIG.*cos(THETg); Y3D = RBIG.*sin(THETg); Z3D = Zg;
            Z3D = Z3D + botExtent-JOINTCUTOFFBOT;
    
            n1 = size(X3D,2); % Number of nodes circumferentially
            n2 = size(X3D,1); % Number of nodes vertically
            N1(S) = n1; N2(S) = n2;
    
            % Storing reconstruction outputs
            CANS.(sprintf('S%g',S)) = [];
    
            CANS.(sprintf('S%g',S)).THETg = THETg;
            CANS.(sprintf('S%g',S)).Zg = Zg;
            CANS.(sprintf('S%g',S)).RBIG = RBIG;
    
            CANS.(sprintf('S%g',S)).X3D = X3D;
            CANS.(sprintf('S%g',S)).Y3D = Y3D;
            CANS.(sprintf('S%g',S)).Z3D = Z3D;
    
            CANS.(sprintf('S%g',S)).N1 = n1;
            CANS.(sprintf('S%g',S)).N2 = n2;
    
        else
            N1(S) = 0; N2(S) = 0;
    
            CANS.(sprintf('S%g',S)) = [];
    
            CANS.(sprintf('S%g',S)).THETg = [];
            CANS.(sprintf('S%g',S)).Zg = [];
            CANS.(sprintf('S%g',S)).RBIG = [];
    
            CANS.(sprintf('S%g',S)).X3D = [];
            CANS.(sprintf('S%g',S)).Y3D = [];
            CANS.(sprintf('S%g',S)).Z3D = [];
    
            CANS.(sprintf('S%g',S)).N1 = 0;
            CANS.(sprintf('S%g',S)).N2 = 0;
        end
    
        disp(' ')
        
    end
    
    % Plotting all reconstructed cans on one figure
    if plotting
        figure
        hold on
        for S = cloudStrakes
            surf(CANS.(sprintf('S%g',S)).THETg,CANS.(sprintf('S%g',S)).Z3D,CANS.(sprintf('S%g',S)).RBIG)
        end
        grid on
    end
    cacheConfig = stage5Config;
    save(s5Cache,'cloudStrakes','CANS','N1','N2','Z_BOTS','Z_TOPS','R0_BOTS','R0_TOPS','THICKS','stageDiagnostics','cacheConfig','-v7.3')
end
if isempty(cloudStrakes)
    error('protocolRun:NoReconstructedStrakes', ...
        'No strakes were reconstructed. Stage 5 saved the exclusion reasons; Stage 6 cannot generate a mesh.')
end

% Produce Stage 5 diagnostics immediately, including the measured cloud.
if generateDiagnostics
    try
        generateProtocolDiagnostics(tower,minZ,maxZ,forceDiagnostics,5)
    catch diagnosticError
        warning('protocolRun:Stage5DiagnosticGenerationFailed', ...
            'Stage 5 outputs are saved, but diagnostic generation failed:\n%s', ...
            getReport(diagnosticError,'extended','hyperlinks','off'))
    end
end

%% STAGE 6: PROJECTION TO ARBITRARY MESHES
% If a later start stage is requested and the relevant .mat file exists, load the .mat file, else run the stage
s6Cache = fullfile(inputDir,['S6_Mesh_',tower,'.mat']);
stage6Config = struct('stage',6,'implementationVersion',implementationVersion, ...
    'upstream',stage5Config);
if startFromStage > 6 && isCompatibleCache(s6Cache,stage6Config)
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Loading')
    loaded = load(s6Cache,'cloudStrakes','X_MESH','Y_MESH','Z_MESH','THET_MESH','R_MESH','N1_MESH','N2_MESH','Z_BOTS','Z_TOPS','R0_BOTS','R0_TOPS','THICKS');
    cloudStrakes = loaded.cloudStrakes; X_MESH = loaded.X_MESH;
    Y_MESH = loaded.Y_MESH; Z_MESH = loaded.Z_MESH;
    THET_MESH = loaded.THET_MESH; R_MESH = loaded.R_MESH;
    N1_MESH = loaded.N1_MESH; N2_MESH = loaded.N2_MESH;
    Z_BOTS = loaded.Z_BOTS; Z_TOPS = loaded.Z_TOPS;
    R0_BOTS = loaded.R0_BOTS; R0_TOPS = loaded.R0_TOPS; THICKS = loaded.THICKS;
else
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Start')
    [X_MESH, Y_MESH, Z_MESH, THET_MESH, R_MESH, N1_MESH, N2_MESH] = S6_ProjectionToArbitraryMeshes(CANS, cloudStrakes, Z_BOTS, Z_TOPS, R0_BOTS, R0_TOPS, THICKS, plotting);
    cacheConfig = stage6Config;
    save(s6Cache,'cloudStrakes','X_MESH','Y_MESH','Z_MESH','THET_MESH','R_MESH','N1_MESH','N2_MESH','Z_BOTS','Z_TOPS','R0_BOTS','R0_TOPS','THICKS','cacheConfig','-v7.3')
    disp('STAGE 6: PROJECTION TO ARBITRARY MESHES - Complete')
end

%% STANDARD DIAGNOSTIC FIGURES
% Stage 5 figures were generated above. Finish the remaining diagnostics now
% that the Stage 6 mesh cache is available.
if generateDiagnostics
    disp('STANDARD DIAGNOSTIC FIGURES - Start')
    try
        generateProtocolDiagnostics(tower,minZ,maxZ,forceDiagnostics,[1:4 6])
        disp('STANDARD DIAGNOSTIC FIGURES - Complete')
    catch diagnosticError
        warning('protocolRun:DiagnosticGenerationFailed', ...
            'Protocol outputs are complete, but diagnostic generation failed:\n%s', ...
            getReport(diagnosticError,'extended','hyperlinks','off'))
    end
end
