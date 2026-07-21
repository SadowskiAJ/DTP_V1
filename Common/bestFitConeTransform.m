function [X, Y, Z, T, R] = bestFitConeTransform(X, Y, Z, minZ, D, plotting)

% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.

% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git

% Function applies best fit cone transformation parameters.

if nargin < 6
    plotting = false;
end

X = reshape(X, length(X), 1);
Y = reshape(Y, length(X), 1);
Z = reshape(Z, length(X), 1);

Z = Z - minZ;
coords = [X-D(1) Y-D(2) Z]; % Shifting origin
Rx = [1 0 0; 0 cos(D(3)) -sin(D(3)); 0 sin(D(3)) cos(D(3))]; % Rotation about x axis
Ry = [cos(D(4)) 0 sin(D(4)); 0 1 0; -sin(D(4)) 0 cos(D(4))]; % Rotation about Y axis

coords = coords*Rx*Ry; % Rotating coordinates

X = reshape(coords(:,1), 1, length(X));
Y = reshape(coords(:,2), 1, length(X));
Z = reshape(coords(:,3), 1, length(X));

Z = Z + minZ;

if plotting
    figure
    hold on
    scatter3(X, Y, Z, 0.1, '.k')
    xlabel('$x$ [m]', 'Interpreter','latex')
    ylabel('$y$ [m]', 'Interpreter','latex')
    zlabel('$z$ [m]', 'Interpreter','latex')
    title('Point cloud after best fit correction', 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    xline(0)
    yline(0)
    axis equal tight
    grid on
end

% Converting to cylindrical coordinates
T = atan2(Y, X);
R = sqrt(Y.^2 + X.^2);
T(T < 0) = T(T < 0) + 2*pi;

end
