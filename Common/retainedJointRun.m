function retained = retainedJointRun(joints,trimFailedEndJoints,retainJointRuns,isFlangeJoint,nominalZ)
% Keep separate usable joint runs, or apply the original end-only trimming.
validateattributes(trimFailedEndJoints,{'logical'},{'scalar'})
if nargin < 3; retainJointRuns = false; end
if nargin < 4; isFlangeJoint = false(size(joints)); end
if nargin < 5; nominalZ = joints; end
validateattributes(retainJointRuns,{'logical'},{'scalar'})
if retainJointRuns
    retained = ~isnan(joints);
    % A flange fit needs its actual neighbouring welds, not the next found peaks.
    for j = find(isFlangeJoint(:))'
        if j == 1 || j == numel(joints) || ~retained(j-1) || ~retained(j+1) || ...
                ~flangeWeldsSupported(joints(j-1:j+1),nominalZ(j-1:j+1))
            retained(j) = false;
        end
    end
    % An isolated joint cannot bound a strake.
    adjacent = retained(1:end-1) & retained(2:end);
    retained = reshape([adjacent(:);false] | [false;adjacent(:)],size(joints));
    return
end
retained = true(size(joints));
if trimFailedEndJoints
    present = find(~isnan(joints));
    retained(:) = false;
    if ~isempty(present)
        retained(present(1):present(end)) = true;
    end
end
end
