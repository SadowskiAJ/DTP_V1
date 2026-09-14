function [joints,matches,selected] = selectNominalJointPeaks(peakZ,peakStrength,nominalZ,radius)
% Select the strongest existing candidate near each nominal junction.
validateattributes(radius,{'numeric'},{'real','scalar','finite','positive'})
nominalZ = nominalZ(:); peakZ = peakZ(:); peakStrength = peakStrength(:);
if isempty(nominalZ) || any(~isfinite(nominalZ)) || any(diff(nominalZ) <= 0) || ...
        numel(peakZ) ~= numel(peakStrength) || any(~isfinite(peakZ))
    error('selectNominalJointPeaks:InvalidInput','Nominal elevations must increase; candidate elevations and strengths must have equal lengths.')
end
midpoints = (nominalZ(1:end-1)+nominalZ(2:end))/2;
lower = max(nominalZ-radius,[-Inf;midpoints]);
upper = min(nominalZ+radius,[midpoints;Inf]);
joints = nan(size(nominalZ)); selected = zeros(size(nominalZ));
strength = nan(size(nominalZ)); candidateCount = zeros(size(nominalZ));
for j = 1:numel(nominalZ)
    candidates = find(peakZ > lower(j) & peakZ < upper(j) & isfinite(peakStrength));
    candidateCount(j) = numel(candidates);
    if isempty(candidates); continue; end
    [strength(j),best] = max(peakStrength(candidates));
    selected(j) = candidates(best);
    joints(j) = peakZ(selected(j));
end
matches = table(nominalZ,joints,joints-nominalZ,lower,upper,candidateCount,strength, ...
    'VariableNames',{'Nominal_m','Detected_m','Offset_m','SearchLower_m','SearchUpper_m','Candidates','PeakStrength_m'});
joints = joints';
end
