function [binParents, binCorners, pointBins] = createOctTree(points, npointsbin, initialBounds)

    % Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
    % of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
    % Energy.
    
    % Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
    
    % Function to be called by user to create the octree
    
    if size(points,2) ~= 3 || isempty(points) || any(~isfinite(points),'all')
        error('createOctTree:InvalidPoints','points must be a finite nonempty N-by-3 array.')
    end
    if ~(isscalar(npointsbin) && npointsbin >= 1 && npointsbin == floor(npointsbin))
        error('createOctTree:InvalidLeafCapacity','npointsbin must be a positive integer.')
    end
    if numel(initialBounds) ~= 6 || any(~isfinite(initialBounds)) || ...
            any(initialBounds(4:6) <= initialBounds(1:3))
        error('createOctTree:InvalidBounds','initialBounds must define a finite positive-volume box.')
    end
    if any(points < initialBounds(1:3) | points > initialBounds(4:6),'all')
        error('createOctTree:PointsOutsideBounds','All points must lie inside initialBounds.')
    end

    pointBins = ones(size(points,1),1); % Stores which bin each point belongs to
    binParents = [0]; % Storing where each point came from
    binCorners = initialBounds;
    binNo = 1;
    [binParents, binCorners, pointBins] = OctTreeRecursive(binNo, binParents, binCorners, pointBins, points, npointsbin, 0);
end

% Funtcion called by createOctTree to recursively create the tree
function [binParents, binCorners, pointBins] = OctTreeRecursive(binNo, binParents, binCorners, pointBins, points, npointsbin, depth)

if depth >= 64
    warning('createOctTree:MaximumDepth','Maximum octree depth reached in bin %g.',binNo)
    return
end

binXl = (binCorners(binNo,4) - binCorners(binNo,1))/2;
binYl = (binCorners(binNo,5) - binCorners(binNo,2))/2;
binZl = (binCorners(binNo,6) - binCorners(binNo,3))/2;
newBinNo = size(binCorners,1); % Updating the newBinNo to the last bin created

for j = 1:8
    if j==1
        binCorner = [binCorners(binNo,1) binCorners(binNo,2) binCorners(binNo,3)+binZl binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,6)];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==2
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2) binCorners(binNo,3)+binZl binCorners(binNo,4) binCorners(binNo,2)+binYl binCorners(binNo,6)];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==3
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2) binCorners(binNo,3) binCorners(binNo,4) binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==4
        binCorner = [binCorners(binNo,1) binCorners(binNo,2) binCorners(binNo,3) binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>=binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==5
        binCorner = [binCorners(binNo,1) binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl binCorners(binNo,1)+binXl binCorners(binNo,5) binCorners(binNo,6)];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==6
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3)+binZl binCorners(binNo,4) binCorners(binNo,5) binCorners(binNo,6)];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<=binCorner(6);
    elseif j==7
        binCorner = [binCorners(binNo,1)+binXl binCorners(binNo,2)+binYl binCorners(binNo,3) binCorners(binNo,4) binCorners(binNo,5) binCorners(binNo,3)+binZl];
        sel = points(:,1)>binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    elseif j==8
        binCorner = [binCorners(binNo,1) binCorners(binNo,2)+binYl binCorners(binNo,3) binCorners(binNo,1)+binXl binCorners(binNo,5) binCorners(binNo,3)+binZl];
        sel = points(:,1)>=binCorner(1) & points(:,1)<=binCorner(4) & points(:,2)>binCorner(2) & points(:,2)<=binCorner(5) & points(:,3)>=binCorner(3) & points(:,3)<binCorner(6);
    end


        newBinNo = newBinNo + 1; % Increasing the bin number if there is anything in this bin
        binParents(newBinNo,1) = binNo; % Storing where the new bins came from
        binCorners(newBinNo,1:6) = binCorner; % Storing the corners of each bin


        pointBins(sel) = newBinNo; % Updating which subdivision a point lies in

        if nnz(sel) > npointsbin
            selectedPoints = points(sel,:);
            coordinateRange = max(selectedPoints,[],1)-min(selectedPoints,[],1);
            scale = max(1,max(abs(selectedPoints),[],'all'));
            if any(coordinateRange > 64*eps(scale))
                [binParents, binCorners, pointBins] = OctTreeRecursive(newBinNo,binParents, binCorners, pointBins, points, npointsbin, depth+1);
            else
                warning('createOctTree:CoincidentPoints', ...
                    'Bin %g contains more than %g numerically coincident points.',newBinNo,npointsbin)
            end
            newBinNo = size(binCorners,1); % Updating the newBinNo to the last bin created
        end

end
end
