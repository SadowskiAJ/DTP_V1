% Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
% of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
% Energy.
%
% Copyright under a BSD 3-Clause License, see
% https://github.com/SadowskiAJ/DTP_V1.git
%
function points = lineCircleSegmentIntersections(p0, p1, centre, radius)
%LINECIRCLESEGMENTINTERSECTIONS Intersections of a segment and a circle.

p0 = reshape(p0, 1, 2);
p1 = reshape(p1, 1, 2);
centre = reshape(centre, 1, 2);
direction = p1-p0;
a = dot(direction, direction);
if a <= eps(max(1, dot(p0, p0)+dot(p1, p1))) || radius < 0
    points = zeros(0,2);
    return
end

offset = p0-centre;
b = 2*dot(offset, direction);
c = dot(offset, offset)-radius^2;
discriminant = b^2-4*a*c;
discTol = 64*eps(max([1, abs(b^2), abs(4*a*c)]));
if discriminant < -discTol
    points = zeros(0,2);
    return
elseif abs(discriminant) <= discTol
    parameters = -b/(2*a);
else
    rootDisc = sqrt(max(discriminant,0));
    parameters = [(-b-rootDisc)/(2*a), (-b+rootDisc)/(2*a)];
end

parameterTol = 64*eps;
parameters = parameters(parameters >= -parameterTol & parameters <= 1+parameterTol);
parameters = min(max(parameters,0),1);
parameters = unique(parameters, 'stable');
points = p0 + parameters(:).*direction;
end
