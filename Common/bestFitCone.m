% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function D = bestFitCone(X, Y, Z, minZ, tower, nRandInds)
%BESTFITCONE Return the parameters of a best-fitting truncated cone.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repoRoot, ['Input_', tower], ['S0_Data_', tower]))

X = X(:);
Y = Y(:);
Z = Z(:) - minZ;
if numel(X) ~= numel(Y) || numel(X) ~= numel(Z) || isempty(X)
    error('bestFitCone:InvalidInput', ...
        'X, Y and Z must be nonempty vectors of equal length.')
end

nSample = min(numel(X), max(1, floor(nRandInds)));
if nSample < numel(X)
    sampleIndices = randperm(numel(X), nSample);
else
    sampleIndices = 1:numel(X);
end
xm = [X(sampleIndices), Y(sampleIndices), Z(sampleIndices)];

options = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', 5e5, ...
    'FunctionTolerance', 1e-14, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 5000);

bottomRadius = R0_BOTS(1)/1e3;
topRadius = R0_TOPS(end)/1e3;
b0 = [mean(X(sampleIndices)), mean(Y(sampleIndices)), 0, 0, ...
    bottomRadius, topRadius];
D = lsqnonlin(@(b) coneFit(b, xm), b0, [], [], options);
end

function perpImp = coneFit(b, x)
X = x(:,1);
Y = x(:,2);
Z = x(:,3);
H = max(Z);
if H <= 0
    error('bestFitCone:InvalidHeight', 'The fitted elevation range must be positive.')
end

coords = [X-b(1), Y-b(2), Z];
Rx = [1 0 0; 0 cos(b(3)) -sin(b(3)); 0 sin(b(3)) cos(b(3))];
Ry = [cos(b(4)) 0 sin(b(4)); 0 1 0; -sin(b(4)) 0 cos(b(4))];
coords = coords*Rx*Ry;

targetRadius = (b(5)-b(6))/H*(H-Z)+b(6);
beta = atan2(b(5)-b(6), H);
radialImp = targetRadius - hypot(coords(:,1), coords(:,2));
perpImp = radialImp.*cos(beta);
end
