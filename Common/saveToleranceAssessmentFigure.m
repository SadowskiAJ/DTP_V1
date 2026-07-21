% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function saveToleranceAssessmentFigure(figureHandle,inputDir,stem)
%SAVETOLERANCEASSESSMENTFIGURE Save an editable and publication-quality figure.
%   The .fig and 300 dpi .png files are written to the corresponding
%   Input_<tower> directory so that tolerance-assessment results remain
%   beside the protocol inputs and cached outputs from which they derive.

if ~isfolder(inputDir)
    error('saveToleranceAssessmentFigure:MissingInputDirectory', ...
        'Input directory does not exist: %s',inputDir)
end
if ~isgraphics(figureHandle,'figure')
    error('saveToleranceAssessmentFigure:InvalidFigure', ...
        'figureHandle must refer to a valid MATLAB figure.')
end

set(figureHandle,'Color','w','InvertHardcopy','off')

axesHandles = findall(figureHandle,'-isa','matlab.graphics.axis.Axes');
for axesIndex = 1:numel(axesHandles)
    axesHandle = axesHandles(axesIndex);
    if ~isappdata(axesHandle,'PreserveToleranceAxesBackground')
        axesHandle.Color = 'w';
    end
    set(axesHandle,'XColor','k','YColor','k','ZColor','k', ...
        'GridColor',[0.75 0.75 0.75], ...
        'MinorGridColor',[0.85 0.85 0.85], ...
        'TickLabelInterpreter','latex','FontSize',12,'LineWidth',0.8)
    axesHandle.Title.Color = 'k';
    axesHandle.Title.Interpreter = 'latex';
    axesHandle.XLabel.Color = 'k';
    axesHandle.XLabel.Interpreter = 'latex';
    axesHandle.YLabel.Color = 'k';
    axesHandle.YLabel.Interpreter = 'latex';
    axesHandle.ZLabel.Color = 'k';
    axesHandle.ZLabel.Interpreter = 'latex';
end

colourBars = findall(figureHandle,'-isa', ...
    'matlab.graphics.illustration.ColorBar');
for colourBarIndex = 1:numel(colourBars)
    colourBar = colourBars(colourBarIndex);
    colourBar.Label.Interpreter = 'latex';
    colourBar.TickLabelInterpreter = 'latex';
    colourBar.FontSize = 12;
    colourBar.Color = 'k';
end

legends = findall(figureHandle,'-isa', ...
    'matlab.graphics.illustration.Legend');
for legendIndex = 1:numel(legends)
    set(legends(legendIndex),'Interpreter','latex','FontSize',11, ...
        'Color','w','TextColor','k','EdgeColor',[0.25 0.25 0.25])
end
drawnow

figFile = fullfile(inputDir,[stem,'.fig']);
pngFile = fullfile(inputDir,[stem,'.png']);
savefig(figureHandle,figFile)
exportgraphics(figureHandle,pngFile, ...
    'Resolution',300,'BackgroundColor','white')

fprintf('Saved tolerance-assessment figures:\n  %s\n  %s\n',figFile,pngFile)
end
