% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function figureHandle = plotToleranceHeatmap(theta,z,radius,values, ...
    plotTitle,colourBarLabel)
%PLOTTOLERANCEHEATMAP Plot a tolerance field over the developed shell.
%   The default view is a clear theta-elevation heat map. The underlying
%   surface and editable .fig remain three-dimensional and rotatable.

finiteValues = values(isfinite(values));
if isempty(finiteValues)
    colourLimit = 1;
else
    colourLimit = max(finiteValues,[],'all');
    colourLimit = max(colourLimit,eps);
end

figureHandle = figure('Color','w','Position',[100 100 1250 720], ...
    'Name',plotTitle,'NumberTitle','off');
axesHandle = axes(figureHandle);
surf(axesHandle,theta,z/1000,radius,values,'EdgeColor','none');
assessedCoverage = 100*nnz(isfinite(values))/numel(values);
if any(~isfinite(values(:)))
    axesHandle.Color = [0.86 0.86 0.86];
    setappdata(axesHandle,'PreserveToleranceAxesBackground',true)
    subtitle(axesHandle, ...
        sprintf(['Assessed coverage: %.1f\\%%; grey denotes locations ' ...
        'not assessed by the gauge algorithm'],assessedCoverage), ...
        'Interpreter','latex')
end
fprintf('Tolerance-map assessed coverage: %.1f%%.\n',assessedCoverage)
view(axesHandle,2)
axis(axesHandle,'tight');
box(axesHandle,'on');
grid(axesHandle,'on');
xlim(axesHandle,[0 2*pi]);
xticks(axesHandle,[0 pi/2 pi 3*pi/2 2*pi]);
xticklabels(axesHandle,{'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'});
clim(axesHandle,[0 colourLimit]);
colormap(axesHandle,turbo)
xlabel(axesHandle,'Circumferential coordinate $\theta$ [rad]', ...
    'Interpreter','latex')
ylabel(axesHandle,'Elevation $z$ [m]','Interpreter','latex')
title(axesHandle,plotTitle,'Interpreter','latex')
colourBar = colorbar(axesHandle);
colourBar.Label.String = colourBarLabel;
colourBar.Label.Interpreter = 'latex';
end
