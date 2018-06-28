function linked_segments = LinkSegs(Segments)
StartPoints = [Segments(:,1),Segments(:,2)];
EndPoints   = [Segments(:,3),Segments(:,4)];
[ValidPts, Distance] = knnsearch(StartPoints,EndPoints,'K',50);
linked_segments = struct([]);
for i = 1:size(Segments,1)
   tmp = CheckSuccessor(i,ValidPts(i,:)',Distance(i,:)',Segments);
   linked_segments = [linked_segments; tmp];
end
end

function LinkedSegments = CheckSuccessor(Pos,PossibleIdx, Distances, SegmentList)
LinkedSegments = struct('SegNum',Pos,'LSeg',[]);

GoodPts = PossibleIdx(Distances < SegmentList(Pos,6));
ValidPts = [];
for i = 1:size(GoodPts)
    if(mod2pi(SegmentList(GoodPts(i),5) - SegmentList(Pos,5)) > 0)
        ValidPts = [ValidPts,GoodPts(i)];
    end
end

if(isempty(ValidPts))
    LinkedSegments.SegNum = Pos;
    return;
end

ParentLine = SegmentList(Pos,1:4);
ChildLine = SegmentList(ValidPts,1:4);

IntersectionPts = [];
for j = 1:size(ValidPts,2)
    [tmpX,tmpY]= IntersectionWith(ParentLine,ChildLine(j,:));
    IntersectionPts = [IntersectionPts;tmpX,tmpY];
end

ParentDist = Pt2PtDist(ParentLine(3),ParentLine(4),IntersectionPts(:,1),IntersectionPts(:,2));
ChildDist  = Pt2PtDist(ChildLine(:,1),ChildLine(:,2),IntersectionPts(:,1),IntersectionPts(:,2));

for i = 1:size(ParentDist,2)
    if(max(ParentDist(i),ChildDist(i)) < SegmentList(Pos,6))
        LinkedSegments.LSeg = [LinkedSegments.LSeg, ValidPts(i)];
    end
end

end

function [x, y] = IntersectionWith(ParentLine, ChildLine)

m00 = ParentLine(1) - ParentLine(3);
m01 = -(ChildLine(1) - ChildLine(3));
m10 = ParentLine(2) - ParentLine(4);
m11 = -(ChildLine(2) - ChildLine(4));

det = m00*m11 - m01*m10;

if(abs(det) < 1e-10)
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
end

function distance = Pt2PtDist(P1x,P1y,P2x,P2y)
dx = P1x - P2x; %Change in X
dy = P1y - P2y; %Change in Y
distance = sqrt(dx.^2 + dy.^2); %Find the Euclidean distance 
end