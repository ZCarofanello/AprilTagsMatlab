function linked_segments = LinkSegs(Segments)

StartPoints = [Segments(:,1),Segments(:,2)]; %Get all the first points of line segments
EndPoints   = [Segments(:,3),Segments(:,4)]; %Get all the second points of line segments

ValidPts = knnsearch(StartPoints,EndPoints,'K',10); %Find the nearest 20 start points to each end point
linked_segments = struct([]); %Allocating a struct to add segments to
for i = 1:size(Segments,1)
   tmp = CheckSuccessor(i,ValidPts(i,:)',Segments); %Checks if the segment has any successors
   linked_segments = [linked_segments; tmp];        %Adds it to the list
end
end

function LinkedSegments = CheckSuccessor(Pos,PossibleIdx, SegmentList)
LinkedSegments = struct('SegNum',Pos,'LSeg',[]);

GoodPts = PossibleIdx;

ValidPts = [];
for i = 1:size(GoodPts)
    if(mod2pi(SegmentList(GoodPts(i),5) - SegmentList(Pos,5)) < 0)
        ValidPts = [ValidPts,GoodPts(i)];
    end
end

if(isempty(ValidPts)) %Oops we don't have any successors
    return;
end

ParentLine = [SegmentList(Pos,1),SegmentList(Pos,2); SegmentList(Pos,3), SegmentList(Pos,4)];     %Pair of x/y coordinates of Parent Line
ChildLine =  SegmentList(ValidPts,1:4); %Pair of x/y coordinates of Child Line

IntersectionPts = []; %Empty matrix
POI = [5,50,6,3];

for j = 1:size(ValidPts,2)
    CurrentChild = [ChildLine(j,1),ChildLine(j,2);ChildLine(j,3),ChildLine(j,4)];
    if(ismember(Pos,POI))
    	tmp = DiffIntersect([ParentLine;CurrentChild]); %calculates the intersections
        hold on
        scatter(tmp(1),tmp(2))
        hold off
    else
        tmp = DiffIntersect([ParentLine;CurrentChild]); %calculates the intersections
    end
    IntersectionPts = [IntersectionPts;tmp];   %Adds that point to list of intersections
end

ParentDist = Pt2PtDist(ParentLine(2,1),ParentLine(2,2),IntersectionPts(:,1),IntersectionPts(:,2));
ChildDist  = Pt2PtDist(ChildLine(:,1),ChildLine(:,2),IntersectionPts(:,1),IntersectionPts(:,2));

for i = 1:size(ParentDist,1)
    if(max(ParentDist(i),ChildDist(i)) < SegmentList(Pos,6))
        LinkedSegments.LSeg = [LinkedSegments.LSeg, ValidPts(i)];
    end
end

end

function [x, y] = IntersectionWith(ParentLine, ChildLine, Debug)

m00 = -(ParentLine(3) - ParentLine(1));  % dx of parent line
m01 = -(ChildLine(1) - ChildLine(3)); % dx of child line
m10 = -(ParentLine(4) - ParentLine(2));  %dy of parent line
m11 = -(ChildLine(2) - ChildLine(4)); %dy of child line

det = m00*m11 - m01*m10; %find the determinant

if(abs(det) < 1e-10) %Oops they're parallel
    x = NaN;
    y = NaN;
    return;
end

i00 = m11/det;
i01 = m01/det;

b00 = ChildLine(1) - ParentLine(1);
b10 = ChildLine(2) - ParentLine(2);

x00 = i00*b00 + i01*b10;

x = m00*x00+ParentLine(1);
y = m10*x00+ParentLine(2);

if(Debug == 1)
    hold on;
    scatter(x,y)
    hold off;
end
end


function point = DiffIntersect(lines)
% calculate intersection point of two 2d lines specified with 2 points each
% (X1, Y1; X2, Y2; X3, Y3; X4, Y4), while 1&2 and 3&4 specify a line.
% Gives back NaN or Inf/-Inf if lines are parallel (= when denominator = 0)
% see http://en.wikipedia.org/wiki/Line-line_intersection
    x = lines(:,1);
    y = lines(:,2);
    % Calculation
    denominator = (x(1)-x(2))*(y(3)-y(4))-(y(1)-y(2))*(x(3)-x(4));
    point = [((x(1)*y(2)-y(1)*x(2))*(x(3)-x(4))-(x(1)-x(2))*(x(3)*y(4)-y(3)*x(4)))/denominator ...
        ,((x(1)*y(2)-y(1)*x(2))*(y(3)-y(4))-(y(1)-y(2))*(x(3)*y(4)-y(3)*x(4)))/denominator];
end

function distance = Pt2PtDist(P1x,P1y,P2x,P2y)
dx = P1x - P2x; %Change in X
dy = P1y - P2y; %Change in Y
distance = sqrt(dx.^2 + dy.^2); %Find the Euclidean distance 
end