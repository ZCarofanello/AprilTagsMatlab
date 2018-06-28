clc;
clear all;
close all;

bike = '../pics/bikesgray.jpg';
tag = '../pics/test_tag.png';
ref = '../pics/tag_middle.png';
real = '../pics/real_life_tag.png';


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
tStep2 = toc(tStart)

%Stage 3 + 4: Edge Extraction / Clustering
image_clusters = CalcEdges(gm,gd,50);
tStep3_4 = toc(tStart) - tStep2
%Stage 5: Segmentation (Need to add the 
MinCluster = 4;
FoundSegs   = Segmenter(image_clusters,gd,image_gray);
tStep5 = toc(tStart) - tStep3_4

%Stage 6: Chain Segments
linked_segments = LinkSegs(FoundSegs);

tElapsed = toc(tStart)

function linked_segments = LinkSegs(Segments)

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

x = dx*x00+ParentLine(1);
y = dy*x00+ParentLine(2);
end

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