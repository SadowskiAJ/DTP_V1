% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function [centre, radius] = bestFitCircle2D(X, Y)
%BESTFITCIRCLE2D Geometric least-squares fit of a circle to planar points.

X = X(:);
Y = Y(:);
if numel(X) ~= numel(Y) || numel(X) < 3 || any(~isfinite([X; Y]))
    error('bestFitCircle2D:InvalidInput', ...
        'X and Y must contain at least three finite coordinate pairs.')
end

x0 = mean(X);
y0 = mean(Y);
r0 = mean(hypot(X-x0, Y-y0));
options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
solution = lsqnonlin(@(b) hypot(X-b(1), Y-b(2))-b(3), ...
    [x0, y0, r0], [-inf, -inf, 0], [], options);
centre = solution(1:2);
radius = solution(3);
end
