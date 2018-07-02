clc;
clear all;
close all;

bike  = '../pics/bikesgray.jpg';
tag   = '../pics/test_tag.png';
ref   = '../pics/tag_middle.png';
real  = '../pics/real_life_tag.png';
real2 = '../pics/real_life_tag2.jpg';

Debug_Gradient = 0;
tStart = tic;
image = imread(tag);
figure('Name','Original Image');
imshow(image);
title('Original Image');

%Preprocessing to Grayscale
if(ndims(image) > 2)
    image_gray = rgb2gray(image);
else
    image_gray = image;
end

%image_gray = NormalizeVals(image_gray);
figure('Name','Preprocessing: Grayscale');
imshow(image_gray);
title('Preprocessing: Grayscale');

%Stage 1: Gaussian Blurring 
image_blurred = imgaussfilt(image_gray, 2.4);
%image_blurred = imgaussfilt(image_blurred, 2.4);
figure('Name','Stage 1:Gaussian Blurring');
imshow(image_blurred);
title('Stage 1:Gaussian Blurring');

%Stage 2: Calculating Gradients (Without toolbox)
G = fspecial('gaussian',3,0.5); %Generate Gausian Filter
[dx, dy] = gradient(G);         %Calc Gradient of Gausian
Ix = conv2(image_blurred,dx,'same');
Iy = conv2(image_blurred,dy,'same');

if(Debug_Gradient == 1)
    Ixn = NormalizeVals(Ix);
    Iyn = NormalizeVals(Iy);
    figure('Name','Stage 2a(Debug): Gradient Magnitue (x direction)');
    imshow(Ixn);
    title('Stage 2a: Gradient Magnitue (x direction)');
    figure('Name','Stage 2a(Debug): Gradient Magnitue (y direction)');
    imshow(Iyn);
    title('Stage 2a: Gradient Magnitue (y direction)');
end

gm = single(sqrt(Ix.^2 + Iy.^2));   %Magnitude
gd = single(atan2(Iy,Ix));          %Direction

%[gm, gd]  = imgradient(image_blurred);


figure('Name','Stage 2a: Gradient Magnitue');
imshow(NormalizeVals(gm));
title('Stage 2a: Gradient Magnitue');
figure('Name','Stage 2b: Gradient Direction');
imshow(NormalizeVals(gd));
title('Stage 2b: Gradient Direction');
tStep2 = toc(tStart)

%Stage 3 + 4: Edge Extraction / Clustering
image_clusters = CalcEdges(gm,gd,5);
tStep3_4 = toc(tStart) - tStep2

Cluster_Num = unique(image_clusters(:,4)); %Gets each unique cluster
current_num = 1; %holds the offset of the where we're grabbing clusters

    figure;
    imshow(image_gray);
    hold on;
    for i = 1:size(Cluster_Num)
        num_of_pts = size(find(image_clusters(:,4) == Cluster_Num(i)),1);
        temp = image_clusters(current_num:num_of_pts+current_num - 1,:);
        plot(temp(:,2),temp(:,1),'*');
        current_num = current_num + num_of_pts; %Add to the offset
    end


%Stage 5: Segmentation 
MinCluster = 4;
FoundSegs   = Segmenter(image_clusters,gd,image_gray);
tStep5 = toc(tStart) - tStep3_4

%Stage 6: Chain Segments
linked_segments = LinkSegs(FoundSegs);
tStep6 = toc(tStart) - tStep5

%Stage 7: Find Quads
quads = QuadDetection(linked_segments,FoundSegs);

tElapsed = toc(tStart)

function [x,y] = OpticalCenter(height,width)
x = round(width/2);
y = round(width/2);
end

function output = NormalizeVals(input,Max,Min)
    switch nargin
        case 1
            output = (input-min(input(:)))./(max(input(:))-min(input(:)));
        otherwise
            output = (input-Min)./(Max-Min);
    end
end

function quads = QuadDetection(LinkSegs,FoundSegs)
for i = 1:length(LinkSegs)
    if(isempty(LinkSegs(i).LSeg))
        continue;
    else
        quads = QuadSearch(LinkSegs(i).SegNum,LinkSegs(i),1,FoundSegs,LinkSegs);
    end
end
end

function Quad = QuadSearch(Path, ParentSeg, Depth, SegmentList,LinkedSegs)
Path(Depth) = ParentSeg.SegNum;
if(Depth > 4)
    Quad = [];
    return;
end

if(Depth == 4)
    disp('We are at depth 4!');
    if(Path(1) == Path(4))
        disp('The first segment is the same as the last one!');
        %Intersections = FindIntersections(Path,SegmentList);
    end
    Quad = [];
    return;
end

FirstAngle = SegmentList(Path(1),5);
for i = 1:length(ParentSeg.LSeg)
    ChildSegment = SegmentList(ParentSeg.LSeg(i),:);
    if(ChildSegment(5) > FirstAngle)
        continue;
    end
    Path(Depth+1) = ParentSeg.LSeg(i);
    Quad = QuadSearch(Path,LinkedSegs(ParentSeg.LSeg(i)),Depth+2,SegmentList,LinkedSegs);
end
    Quad = [];
end

function QuadPts = FindIntersections(Path, SegmentList)
for i = 1:4
    LineA = SegmentList(Path(i),1:4);
    LineB = SegmentList(Path(i+1),1:4);
    Point(i,:) = IntersectionWith(LineA,LineB);
    if(isNaN(Point(i,1)))
        QuadPts = [NaN,Nan;NaN,Nan;NaN,Nan;NaN,NaN];
        return;
    end
end
QuadPts = Point;
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