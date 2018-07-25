function quads = QuadDetection(LinkSegs,FoundSegs)
quads = [];
%Check through all the segments for possible quads
for i = 1:length(LinkSegs)
    if(isempty(LinkSegs(i).LSeg))
        continue;
    else
        quads = [quads; QuadSearch(LinkSegs(i).SegNum,LinkSegs(i)...
            ,1,FoundSegs,LinkSegs,[])];
    end
end
end

function Quad = QuadSearch(Path, ParentSeg, Depth, SegmentList,LinkedSegs,Quad)
Path(Depth) = ParentSeg.SegNum; %sets current path segment
Quad = Quad; %Fix for recursion

%if Depth is 5 then that means we can have a quad with an 'extra' segment
%that is the same as Depth 1
if(Depth == 5)
    %if first point is equal to the last point
    if(Path(1) == Path(5))
        Bad = false;
        
        %find all the intersections between the segments
        Intersections = FindIntersections(Path,SegmentList);
        
        %Make sure that we make a quad and not an hourglass
        Bad = CheckAngles(Intersections);
        
        if(Bad)
            %disp('Failed at CheckAngles');  
            return;
        end
        
        %Make sure it's not too small
        Bad = CheckSize(Intersections);
        
        if(Bad)
            %disp('Failed at CheckSize or CheckAspectRatio');
            return;
        end
        %Get the points that make up the quad
        Qpoints = [Intersections(1,:),Intersections(2,:)...
            ,Intersections(3,:),Intersections(4,:)];
        
        %Add it to list
        Quad = [Quad;Qpoints];
    end
    return;
end

FirstAngle = SegmentList(Path(1),5); %Needed to make sure winding is obeyed

%Go through all children for this segment
for i = 1:length(ParentSeg.LSeg)
    ChildSegment = SegmentList(ParentSeg.LSeg(i),:);
    
    %Make sure it obeys winding
    if(ChildSegment(5) > FirstAngle)
        continue;
    end
    
    %Recurse down and go through it again
    Quad = QuadSearch(Path,LinkedSegs(ParentSeg.LSeg(i)),...
        Depth+1,SegmentList,LinkedSegs,Quad);
end
end

function QuadPts = FindIntersections(Path, SegmentList)
%Point = zeros(1,2);
for i = 1:4
    %Make line A
    LineA = [SegmentList(Path(i),1), SegmentList(Path(i),2)];
    LineA = [LineA; SegmentList(Path(i),3), SegmentList(Path(i),4)];
    
    %Make line B
    LineB = [SegmentList(Path(i+1),1), SegmentList(Path(i+1),2)];
    LineB = [LineB; SegmentList(Path(i+1),3), SegmentList(Path(i+1),4)];

    %Get the intersection between line A and line B
    Point(i,:) = DiffIntersect([LineA;LineB]);
    
    %Oops they're parallel (shouldn't happen because we're looking down
    %the sides of the quad
    if(Point(i,1) == NaN)
        QuadPts = [NaN,NaN;NaN,NaN;NaN,NaN;NaN,NaN];
        return;
    end
end
QuadPts = Point;
end

function Bad = CheckAngles(Intersections)
pts = Intersections;

%Make sure the angles don't make an hourglass shape with the points
t0 = atan2(pts(2,2) - pts(1,2),pts(2,1) - pts(1,1));
t1 = atan2(pts(3,2) - pts(2,2),pts(3,1) - pts(2,1));
t2 = atan2(pts(4,2) - pts(3,2),pts(4,1) - pts(3,1));
t3 = atan2(pts(1,2) - pts(4,2),pts(1,1) - pts(4,1));

%Add all the angles together
theta = mod2pi(t1-t0) + mod2pi(t2-t1) + mod2pi(t3-t2) + mod2pi(t0-t3);

if(theta < -7 || theta > -5)
    Bad = true;
else
    Bad = false;
end

end

function Bad = CheckSize(Intersections)
%Need to Export
MinEdgeLength = 6;
MaxQuadAspectRatio = 32;

pts = Intersections;

%find distance between the points that make up the quad
d0 = Pt2PtDist(pts(1,1),pts(1,2),pts(2,1),pts(2,2));
d1 = Pt2PtDist(pts(2,1),pts(2,2),pts(3,1),pts(3,2));
d2 = Pt2PtDist(pts(3,1),pts(3,2),pts(4,1),pts(4,2));
d3 = Pt2PtDist(pts(4,1),pts(4,2),pts(1,1),pts(1,2));
d4 = Pt2PtDist(pts(1,1),pts(1,2),pts(3,1),pts(3,2));
d5 = Pt2PtDist(pts(2,1),pts(2,2),pts(4,1),pts(4,2));

%Check if all the distances are above the minimum edge length
V1 = d0 < MinEdgeLength;
V2 = d1 < MinEdgeLength;
V3 = d2 < MinEdgeLength;
V4 = d3 < MinEdgeLength;
V5 = d4 < MinEdgeLength;
V6 = d5 < MinEdgeLength;

if(V1 || V2 || V3 || V4 || V5 || V6)
    Bad = true;
else
    %Check Aspect Ratio
    dmax = max(max(d0,d1), max(d2,d3));
    dmin = min(min(d0,d1), min(d2,d3));
    if( dmax > dmin*MaxQuadAspectRatio)
        Bad = true;
    else
        Bad = false;
    end
end

end

function distance = Pt2PtDist(P1x,P1y,P2x,P2y)
dx = P1x - P2x; %Change in X
dy = P1y - P2y; %Change in Y
distance = sqrt(dx.^2 + dy.^2); %Find the Euclidean distance 
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