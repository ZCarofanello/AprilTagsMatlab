clc;
clear all;
close all;

bike = '../pics/bikesgray.jpg';
tag = '../pics/test_tag.png';
ref = '../pics/tag_middle.png';

Debug_Gradient = 0;

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
image_blurred = imgaussfilt(image_gray, 0.8);
figure('Name','Stage 1:Gaussian Blurring');
imshow(image_blurred);
title('Stage 1:Gaussian Blurring');

%Stage 2: Calculating Gradients (Without toolbox)
G = fspecial('gaussian',3,0.5); %Generate Gausian Filter
[dx, dy] = gradient(G);         %Calc Gradient of Gausian
Ix = conv2(image_gray,dx,'same');
Iy = conv2(image_gray,dy,'same');

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

gm = sqrt(Ix.^2 + Iy.^2);   %Magnitude
gd = atan2(Iy,Ix);          %Direction
figure('Name','Stage 2a: Gradient Magnitue');
imshow(NormalizeVals(gm));
title('Stage 2a: Gradient Magnitue');
figure('Name','Stage 2b: Gradient Direction');
imshow(NormalizeVals(gd));
title('Stage 2b: Gradient Direction');



%houghTest(image_gray);

%Stage 3 + 4: Edge Extraction / Clustering
image_clusters = EdgeFunction(gm,gd,50);

%Stage 5: Segmentation (Need to add the 
MinCluster = 4;
Segments   = Segmenter(image_clusters,gd,image_gray);

function output = NormalizeVals(input,Max,Min)
    switch nargin
        case 1
            output = (input-min(input(:)))./(max(input(:))-min(input(:)));
        otherwise
            output = (input-Min)./(Max-Min);
    end
end

function FoundEdges = EdgeFunction(Magnitude, Direction, MagThr)
Magnitude(Magnitude <= MagThr) = 0;%Makes sure all edges are above threshold
figure;
imshow(Magnitude);
FoundEdges = CalcEdges(Magnitude, Direction);
end

function lines = Segmenter(image_clusters,Theta,image_grey)
MinDist = 4;
Cluster_Num = unique(image_clusters(:,4)); %Gets each unique cluster

segments = []; %Array for holding segments
current_num = 1; %holds the offset of the where we're grabbing clusters
figure;
imshow(image_grey);
hold on;
for i = 1:size(Cluster_Num)
    num_of_pts = size(find(image_clusters(:,4) == Cluster_Num(i)),1);
    
    temp = image_clusters(current_num:num_of_pts+current_num - 1,:); %Get all the points in that cluster
    
    LineTemp = Line2D(temp(:,1),temp(:,2),temp(:,3));                %Find the line created by those points
    
    if(MinDist < Pt2PtDist(LineTemp(1,1),LineTemp(1,2),LineTemp(1,3),LineTemp(1,4)))
        segments = [segments;LineTemp]; %Add to the good segments
        plot([LineTemp(1),LineTemp(3)],[LineTemp(2),LineTemp(4)],'-r'); %plot the segment
    end
    
    current_num = current_num + num_of_pts; %Add to the offset
end
hold off;

lines = segments; %Export those segments

end

function line = Line2D(x,y,w)
weightedX = x .* w; %Weights all the x components
weightedY = y .* w; %Weights all the y components

mX = sum(weightedX); %Sums all the weighted x components
mY = sum(weightedY); %Sums all the weighted y components


mXX = sum(weightedX .* x); %Weighted sum of x squares
mYY = sum(weightedY .* y); %Weighted sum of y squares
mXY = sum(weightedY .* x); %Weighted sum of xy product
n = sum(w);                %Sum of weights

Ex  = mX/n;
Ey  = mY/n;
Cxx = (mXX/n) - Ex^2; 
Cyy = (mYY/n) - Ey^2;
Cxy = (mXY/n) - (Ex*Ey);

phi = 0.5*atan2(-2*Cxy,(Cyy-Cxx)); %Uses SVD to find direction

dx = -sin(phi); %Change in x
dy =  cos(phi); %Change in y
xp = round(Ex); %Rounded X point on line
yp = round(Ey); %Rounded Y point on line

[xp,yp,dx,dy] = normalizeP(dx,dy,xp,yp); %Normalize the point

line_coord = GetLineCoord(x,y,dx,dy);    %Get each coordinate on the line from our data set

maxcoord = max(line_coord); %Find the end of the line
mincoord = min(line_coord); %Find the beginning of the line

[max_x, max_y] = GetPtCoord(xp,yp,dx,dy,maxcoord); %Find where the beginning is on the line
[min_x, min_y] = GetPtCoord(xp,yp,dx,dy,mincoord); %Find where the end is on the line

line = [min_x,min_y,max_x,max_y,0,0];
end

function [y,x] = GetPtCoord(xp,yp,dx,dy,coord)
x = round(xp + coord*dx);
y = round(yp + coord*dy);
end

function coord = GetLineCoord(x,y,dx,dy)
coord = x*dx + y*dy;
end

function [xn,yn,dx,dy] = normalizeP(dx,dy,x,y)
[dx,dy] = normalizeSlope(dx,dy);
dotprod = -dy*x + dx*y;

xn = dotprod * -dy;
yn = dotprod * dx;
end

function [dx,dy] = normalizeSlope(dx,dy)
mag = sqrt(dx^2+dy^2);
dx = dx / mag;
dy = dy / mag;
end

function distance = Pt2PtDist(P1x,P1y,P2x,P2y)
dx = P1x - P2x;
dy = P1y - P2y;
distance = sqrt(dx^2 + dy^2);
end
