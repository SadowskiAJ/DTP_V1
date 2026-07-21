% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [X_MESH, Y_MESH, Z_MESH, THET_MESH, R_MESH, N1_MESH, N2_MESH] = S6_ProjectionToArbitraryMeshes(CANS, cloudStrakes, Z_BOTS, Z_TOPS, R0_BOTS, R0_TOPS, THICKS, plotting)
%S6_PROJECTIONTOARBITRARYMESHES Generate a Koiter-ellipse-compliant mesh.

if nargin < 8
    plotting = false;
end
cloudStrakes = cloudStrakes(:)';
if isempty(cloudStrakes) || any(cloudStrakes ~= floor(cloudStrakes)) || ...
        any(cloudStrakes < 1 | cloudStrakes > numel(Z_BOTS)) || ...
        numel(unique(cloudStrakes)) ~= numel(cloudStrakes)
    error('S6_ProjectionToArbitraryMeshes:InvalidStrakes', ...
        'cloudStrakes must contain unique valid integer strake indices.')
end

N1_MESH = zeros(size(Z_BOTS));
N2_MESH = zeros(size(Z_BOTS));
meridionalNodes = zeros(size(Z_BOTS));
circumferentialElements = 0;
ellipseFactor = 1.0;
for S = cloudStrakes
    rb = R0_BOTS(S);
    rt = R0_TOPS(S);
    height = Z_TOPS(S)-Z_BOTS(S);
    thickness = THICKS(S);
    beta = atan(abs(rb-rt)/height);
    averageRadius = (rb+rt)/(2*cos(beta));
    meridionalLength = height/cos(beta);
    circumferentialElements = max(circumferentialElements, ...
        18.18*ellipseFactor*cos(beta)*sqrt(averageRadius/thickness));
    meridionalElements = ceil(5.79*ellipseFactor*meridionalLength/ ...
        sqrt(averageRadius*thickness));
    meridionalNodes(S) = max(meridionalElements,2)+1;
end
circumferentialNodes = max(3,ceil(circumferentialElements));

totalRows = sum(meridionalNodes(cloudStrakes));
X_MESH = zeros(totalRows,circumferentialNodes);
Y_MESH = zeros(totalRows,circumferentialNodes);
Z_MESH = zeros(totalRows,circumferentialNodes);
THET_MESH = zeros(totalRows,circumferentialNodes);
R_MESH = zeros(totalRows,circumferentialNodes);

firstRow = 1;
for S = cloudStrakes
    targetRows = meridionalNodes(S);
    N1_MESH(S) = circumferentialNodes;
    N2_MESH(S) = targetRows;

    fieldName = sprintf('S%g',S);
    X3D = CANS.(fieldName).X3D;
    Y3D = CANS.(fieldName).Y3D;
    Z3D = CANS.(fieldName).Z3D;
    thetaSource = CANS.(fieldName).THETg(1,:);
    zSource = Z3D(:,1);
    sourceRadii = hypot(X3D,Y3D);
    if any(diff(thetaSource) <= 0) || any(diff(zSource) <= 0)
        error('S6_ProjectionToArbitraryMeshes:NonmonotonicGrid', ...
            'The reconstructed grid for strake %g is not strictly monotonic.',S)
    end

    thetaSource = [thetaSource,2*pi];
    sourceRadii = [sourceRadii,sourceRadii(:,1)];
    thetaTarget = (0:circumferentialNodes-1)*(2*pi/circumferentialNodes);
    zTarget = linspace(zSource(1),zSource(end),targetRows);
    [THETgd,Zgd] = meshgrid(thetaTarget,zTarget);

    interpolant = griddedInterpolant({zSource,thetaSource},sourceRadii, ...
        'linear','linear');
    Rgd = interpolant(Zgd,THETgd);
    Xgd = Rgd.*cos(THETgd);
    Ygd = Rgd.*sin(THETgd);

    rows = firstRow:firstRow+targetRows-1;
    X_MESH(rows,:) = Xgd;
    Y_MESH(rows,:) = Ygd;
    Z_MESH(rows,:) = Zgd;
    THET_MESH(rows,:) = THETgd;
    R_MESH(rows,:) = Rgd;
    firstRow = rows(end)+1;

    if plotting
        figure
        tiledlayout(1,2)
        nexttile
        surf(CANS.(fieldName).THETg,Z3D,sourceRadii(:,1:end-1), ...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
        hold on
        surf(THETgd,Zgd,Rgd)
        grid on
        xlabel('$\theta$ [rad]','interpreter','latex')
        ylabel('$z$ [mm]','interpreter','latex')
        zlabel('$\rho$ [mm]','interpreter','latex')
        nexttile
        surf(X3D,Y3D,Z3D,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
        hold on
        surf(Xgd,Ygd,Zgd)
        axis equal tight
        grid on
    end
end

if plotting
    figure
    surf(THET_MESH,Z_MESH,R_MESH)
    grid on
    xlabel('$\theta$ [rad]','interpreter','latex')
    ylabel('$z$ [mm]','interpreter','latex')
    zlabel('$\rho$ [mm]','interpreter','latex')
    figure
    surf(X_MESH,Y_MESH,Z_MESH)
    axis equal tight
    grid on
end
end
