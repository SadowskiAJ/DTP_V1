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

tower = 'ReleaseTower'; % Tower to be processed
startFromStage = 1; % Stage to start protocol from - 1 to 7. (3, 4 and 5 run within the same loop). Use 7 to skip to final protocol output.
randomSeed = 1729; % Reproducible random sampling and global optimisation
plotting = false; % Diagnostic figures are expensive for full-resolution runs
generateDiagnostics = true; % Save standard Stage 1, 2, 5 and 6 FIG/PNG diagnostics
forceDiagnostics = false; % Regenerate diagnostics even when they are newer than their source MAT files
implementationVersion = 2; % Increment when a code change invalidates cached stages
rng(randomSeed,'twister')
repoRoot = fileparts(mfilename('fullpath'));
inputDir = fullfile(repoRoot,['Input_',tower]);

%% Initial set up before registration
minZ = 8; maxZ = 53; % Range of tower being considered between minZ and maxZ
units  = 'm'; % Units to be used in plots

%% Adding folders to path
addpath(genpath(fullfile(repoRoot,'Common')))
addpath(genpath(inputDir))
addpath(genpath(fullfile(repoRoot,'Stage_1')))
addpath(genpath(fullfile(repoRoot,'Stage_2')))
addpath(genpath(fullfile(repoRoot,'Stage_3')))
addpath(genpath(fullfile(repoRoot,'Stage_4')))
addpath(genpath(fullfile(repoRoot,'Stage_5')))
addpath(genpath(fullfile(repoRoot,'Stage_6')))

%% Create .bin files
convertPTS2BIN(tower) % Creates .bin files if no .bin files exist in dataset directory
inputSignature = pointCloudFileSignature(inputDir);

%% STAGE 1: ENHANCED REGISTRATION
% #### Hyperparameters #### %
nperc = 0.01; % Remaining fraction after downsampling
numIter = 5; % Number of iterations
cutoffReg = 0.02; % Cutoff for nearest neighbour [m]
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
% #### Hyperparameters - Global radial profile generation #### %
nperc = 0.1;
expectedBeadWidth = 0.08; % m - Setting expected bead width
numCircWindows = 100; % Number of circumferential windows
subIntervals = 10; % Number of subIntervals. This controls the window size. Must be at least 2.
% iUse = 2; % Computed with smallest 99th percentile point to plane metric. Uncomment to overule and manually choose a different iteration

% #### Hyperparameters - Joint ROI identification (Some hyperparameters carry over from global radial profile) #### %
nJoints = 16; % Number of joints in data

% #### Hyperparameters - Joint segmentation #### %
flangeTransitionJointIndex = [9]; % Indexes of jointCoordinates where flange-flange joint exists
flangePadding = 10; % [mm] Amount to extend region identified as flange-flange joint
weldPadding = 5; % [mm] Amount to extend region identified as welded joint
s2Cache = fullfile(inputDir,['S2_CanSegmentation_',tower,'.mat']);
stage2Config = struct('stage',2,'implementationVersion',implementationVersion, ...
    'upstream',stage1Config,'nperc',nperc,'expectedBeadWidth',expectedBeadWidth, ...
    'numCircWindows',numCircWindows,'subIntervals',subIntervals,'nJoints',nJoints, ...
    'flangeTransitionJointIndex',flangeTransitionJointIndex, ...
    'flangePadding',flangePadding,'weldPadding',weldPadding);

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
    [ZMEDIANS, RMEDIANS] = S2_GlobalRadialProfileGeneration(nperc, minZ, maxZ, bestFitParamsPreReg, regParams, iUse, tower, expectedBeadWidth, numCircWindows, subIntervals, plotting);
    
    % Joint ROI identification
    disp('STAGE 2: Joint ROI identification')
    [jointCoordinates] = S2_JointROIIdentification(ZMEDIANS, RMEDIANS, nJoints, expectedBeadWidth, plotting);

    % Joint segmentation
    disp('STAGE 2: Joint segmentation')
    [JOINTWIDTH, JOINTBOUNDS] = S2_JointSegmentation(ZMEDIANS, RMEDIANS, jointCoordinates, flangeTransitionJointIndex, flangePadding, weldPadding, plotting);
    cacheConfig = stage2Config;
    save(s2Cache,'ZMEDIANS','RMEDIANS','jointCoordinates','JOINTWIDTH','JOINTBOUNDS','cacheConfig','-v7.3');

    disp('STAGE 2: CAN SEGMENTATION - Complete')
end

%% Preparing data for looping through cans
cloudStrakes = 1:15; % Strakes to be meshed from cloud data - Check if within minZ and maxZ
nperc = 0.02; % Fraction of data to use for remaining stages
cutOff = 1; % [mm] Stage 3 outlier cutoff distance
windowWidthSF = 3; % Stage 3 window-width scale factor
p = 1; % Stage 4 inverse-distance power
searchRadiusGridSpacingSF = 4; % Stage 4 initial search-radius scale factor
smoothingFilterStd = 10; % Stage 4 Gaussian standard deviation in grid cells
s5Cache = fullfile(inputDir,['S5_Surface_',tower,'.mat']);
stage5Config = struct('stage',[3,4,5],'implementationVersion',implementationVersion, ...
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
    
    sel = Z_BOTS/1e3 > minZ & Z_BOTS/1e3 < maxZ;
    firstStrakeInData = find(sel,1); % ID of first strake in data i.e bottom of which strake is first found in jointCoordinates
    JOINTBOUNDS = [zeros(firstStrakeInData-1,2); JOINTBOUNDS]; % Adding extra numbers to JOINTBOUNDS to aid indexing and ensure the nth joint corresponds to the nth can if a different region of shell is considered via different minZ and maxZ
    
    % Applying unit conversion of point cloud to mm
    X = X*1000; Y = Y*1000; Z = Z*1000; R = R*1000;
    JOINTBOUNDS = JOINTBOUNDS*1e3;
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
            [THETg, Zg, RBIG] = S4_GriddedOuterShellSurfaceReconstruction(TT_REM, RR_REM, ZZ_REM, botExtent, topExtent, R0_BOTS, R0_TOPS, Z_BOTS, Z_TOPS, THICKS, S, searchRadiusGridSpacingSF, p, smoothingFilterStd, plotting);
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
% Generate diagnostics only after every protocol stage is complete so the
% Stage 6 mesh cache is always available on a fresh run.
if generateDiagnostics
    disp('STANDARD DIAGNOSTIC FIGURES - Start')
    try
        generateProtocolDiagnostics(tower,minZ,maxZ,forceDiagnostics)
        disp('STANDARD DIAGNOSTIC FIGURES - Complete')
    catch diagnosticError
        warning('protocolRun:DiagnosticGenerationFailed', ...
            'Protocol outputs are complete, but diagnostic generation failed:\n%s', ...
            getReport(diagnosticError,'extended','hyperlinks','off'))
    end
end
