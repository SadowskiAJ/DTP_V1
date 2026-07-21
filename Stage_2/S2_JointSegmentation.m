function [JOINTWIDTH, JOINTBOUNDS] = S2_JointSegmentation(ZMEDIANS, RMEDIANS, jointCoordinates, flangeTransitionJointIndex, flangePadding, weldPadding, plotting)

% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Function fits three piece functions to flange-flange transitions and four
% piece functions to welded joints. 

if nargin < 7
    plotting = false;
end
if any(flangeTransitionJointIndex <= 1 | flangeTransitionJointIndex >= numel(jointCoordinates))
    error('S2_JointSegmentation:InvalidFlangeIndex', ...
        'A flange-transition joint must have an identified joint on both sides.')
end

% Additional parameters
distanceFromEndsOfAdjacentFlanges = 0.05; % Distance from start and end of adjacent flanges to define region for fitting the flange-to-flange transition
distanceFromWeldCentre = 0.1; % Points within this distance from the centre of the weld are considered when fitting the weld

%% Looping through joints

JOINTWIDTH = zeros(length(jointCoordinates), 1);
JOINTBOUNDS = zeros(length(jointCoordinates), 2);

% Looping through each identified joint
for j = 1:length(jointCoordinates)

    if ismember(j, flangeTransitionJointIndex)
        % Flange-to-flange region
        sel = ZMEDIANS > jointCoordinates(j-1)+distanceFromEndsOfAdjacentFlanges & ZMEDIANS < jointCoordinates(j+1)-distanceFromEndsOfAdjacentFlanges;
    else
        % Weld region
        sel = ZMEDIANS > jointCoordinates(j)-distanceFromWeldCentre & ZMEDIANS < jointCoordinates(j)+distanceFromWeldCentre;
    end
    Z = ZMEDIANS(sel);
    R = RMEDIANS(sel);
    if numel(Z) < 5
        error('S2_JointSegmentation:InsufficientData', ...
            'Joint %g has too few profile points for segmentation.', j)
    end

    zRange = max(Z) - min(Z);

    if ismember(j, flangeTransitionJointIndex)
        % Flange-to-flange transition fitting
        N = 3;
        [xs,  ~] = pieceWiseGA(N, R, Z, zRange);
    else
        % Weld fitting
        N = 4;
        [xs,  ~] = pieceWiseGA(N, R, Z, zRange);
    end

    % Extracting bounds
    botBound = xs(2);
    topBound = xs(2*N-2);

    % Padding bounds
    if ismember(j, flangeTransitionJointIndex)
        paddedBotBound = botBound-flangePadding/1e3;
        paddedTopBound = topBound+flangePadding/1e3;
    else
        paddedBotBound = botBound-weldPadding/1e3;
        paddedTopBound = topBound+weldPadding/1e3;
    end

    if plotting
        figure('Units', 'pixels','Position', [680 458 560 420])
        hold on
        p0 = scatter(Z,R,100,'.k');
        p1 = plot([min(Z) xs(2:2:2*N-2) max(Z)], [xs(1:2:2*N-1) xs(2*N)], '-x', 'LineWidth',3);
        grid on
        title(['N = ', num2str(N)], 'Interpreter','latex')
        p2 = xline(botBound, 'k');
        xline(topBound, 'k')
        p3 = xline(paddedBotBound, '--k');
        xline(paddedTopBound, '--k');
        xlabel('$z$ [m]', 'Interpreter','latex')
        ylabel('$\rho$ [m]', 'Interpreter','latex')
        legend([p0, p1, p2, p3], {'Fitting data','Fitted function', 'Original bounds', 'Padded bounds'}, 'Interpreter','latex', 'Location','best')
        set(gca,'TickLabelInterpreter','latex')
    end
    disp([num2str(jointCoordinates(j)), ' m, N = ', num2str(N)])

    JOINTWIDTH(j) = topBound-botBound;
    JOINTBOUNDS(j,:) = [paddedBotBound paddedTopBound];
end

%% Plotting fitted joint bounds against global profile
if plotting
    figure;
    plot(ZMEDIANS ,RMEDIANS)
    hold on
    xline(JOINTBOUNDS(:,1), 'b');
    xline(JOINTBOUNDS(:,2), 'r');
    xlabel('$z$ [m]', 'interpreter', 'latex')
    ylabel('$\rho$ [m]', 'interpreter', 'latex')
    title('Identified bounds of joint regions', 'interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
    grid on
    p1 = xline(JOINTBOUNDS(1,1), 'b');
    p2 = xline(JOINTBOUNDS(1,2), 'r');
    legend([p1, p2], {'Start of joint', 'End of joint'}, 'Interpreter','latex')
end


end

%% Functions
function [rsq] = piecewiseR2(b, N, Z, R)

% Function computes R2
targ = zeros(size(Z));
if N == 1
    x1 = min(Z); y1 = b(1);
    x2 = max(Z); y2 = b(2);
    m = (y2-y1)/(x2-x1);
    c = y1 - m*x1;
    targ = m*Z+c;
else
    for NN = 1:N
        if NN == 1
            sel = Z <= b(2);
            x1 = min(Z); y1 = b(1);
            x2 = b(2); y2 = b(3);
            m = (y2-y1)/(x2-x1);
            c = y1 - m*x1;
            targ(sel) = m*Z(sel)+c;
        elseif NN == N
            sel = Z >= b(2*NN-2);
            x1 = b(2*NN-2); y1 = b(2*N-1);
            x2 = max(Z); y2 = b(2*NN);
            m = (y2-y1)/(x2-x1);
            c = y1 - m*x1;
            targ(sel) = m*Z(sel)+c;
        else
            sel = Z >= b(2*NN-2) & Z <= b(2*NN);
            x1 = b(2*NN-2); y1 = b(2*NN-1);
            x2 = b(2*NN); y2 = b(2*NN+1);
            m = (y2-y1)/(x2-x1);
            c = y1 - m*x1;
            targ(sel) = m*Z(sel)+c;
        end
    end
end
yresid = R-targ;

SSresid = sum(yresid.^2);
SStotal = sum((R-mean(R)).^2);
rsq = 1 - SSresid/SStotal;
end

function [xs, fval] = pieceWiseGA(N, Rnew, Znew, zRange)

% Function solves for piecewise function parameters
ObjectiveFunction = @(b)-piecewiseR2(b, N, Znew, Rnew);
ConstraintFunction = @(b)simple_constraint(b, N);

lb = ones(1,2*N)*min(Rnew);
lb(2:2:2*N-2) = min(Znew);

ub = ones(1,2*N)*max(Rnew);
ub(2:2:2*N-2) = max(Znew);

nvars = N*2;

xm = ones(1,2*N)*median(Rnew);
xm(2:2:2*N-2) = min(Znew)+zRange.*(1:N-1)/N;


options = optimoptions("ga");
options.InitialPopulationMatrix = xm;
options.FunctionTolerance = 1e-10;
options.ConstraintTolerance = 1e-12;
options.InitialPopulationRange = [xm-0.01;xm+0.01];


% Solving
[xs,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],lb,ub,ConstraintFunction, options);

fval = -fval;
end

function [c, ceq] = simple_constraint(b,N)

ceq = [];
c = [b(2:2:2*N-4)-b(4:2:2*N-2)]';

end
