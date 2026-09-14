function retained = retainedJointRun(joints,trimFailedEndJoints)
% Trim only missing end junctions. Interior failures remain for validation.
validateattributes(trimFailedEndJoints,{'logical'},{'scalar'})
retained = true(size(joints));
if trimFailedEndJoints
    present = find(~isnan(joints));
    retained(:) = false;
    if ~isempty(present)
        retained(present(1):present(end)) = true;
    end
end
end
