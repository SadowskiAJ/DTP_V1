function validateJointDetection(detected,nominal,tolerance)
% Check the existing detector's output without changing or selecting peaks.
validateattributes(tolerance,{'numeric'},{'real','scalar','finite','positive'})
detected = detected(:); nominal = nominal(:);
if isempty(nominal) || any(~isfinite(nominal)) || any(diff(nominal) <= 0)
    error('validateJointDetection:InvalidNominal','No valid increasing nominal junction list is available for this range.')
end
if numel(detected) ~= numel(nominal)
    error('validateJointDetection:CountMismatch', ...
        'Stage 2 has %g detected junctions but the input geometry specifies %g. Review Stage 2 before reconstruction.', ...
        numel(detected),numel(nominal))
end
if any(isnan(detected))
    disp(table(find(isnan(detected)),nominal(isnan(detected)), ...
        'VariableNames',{'Joint','MissingNominal_m'}))
    error('validateJointDetection:MissingJunction', ...
        '%g nominal junctions have no candidate peak in their search windows. Review the saved Stage 2 diagnostic.',nnz(isnan(detected)))
end
if any(~isfinite(detected)) || any(diff(detected) <= 0)
    error('validateJointDetection:InvalidDetection','Detected junctions must be finite and strictly increasing.')
end
% Do not let the check accept a peak across an adjacent nominal midpoint.
midpoints = (nominal(1:end-1)+nominal(2:end))/2;
lower = max(nominal-tolerance,[-Inf;midpoints]);
upper = min(nominal+tolerance,[midpoints;Inf]);
bad = detected <= lower | detected >= upper;
if any(bad)
    disp(table(find(bad),nominal(bad),detected(bad),detected(bad)-nominal(bad), ...
        'VariableNames',{'Joint','Nominal_m','Detected_m','Offset_m'}))
    error('validateJointDetection:PositionMismatch', ...
        ['%g junctions disagree with nominal geometry (tolerance %.3f m, limited by adjacent midpoints). ', ...
         'Review the profile and input geometry; peaks have not been reassigned.'],nnz(bad),tolerance)
end
end
