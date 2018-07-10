clc;
clear;
close all;

%Read in Reference Values generated by cpp library
RefBw    = csvread('../data/BWPrint.dat');
RefBlur  = csvread('../data/BlurPrint.dat');
RefMag   = csvread('../data/MagPrint.dat');
RefTheta = csvread('../data/ThetaPrint.dat');

bike  = '../pics/bikesgray.jpg';
tag   = '../pics/test_tag.png';
ref   = '../pics/tag_middle.png';
real  = '../pics/real_life_tag.png';
real2 = '../pics/real_life_tag2.jpg';

downscale = 1/255;
Debug_Gradient = 0;
tStart = tic;
image = imread(tag);
figure('Name','Original Image');
imshow(image);
title('Original Image');

%Constants
TagSize = 0.166;
Fx = 600;
Fy = 600;
Px = size(image,2)/2;
Py = size(image,1)/2;

%Preprocessing to Grayscale
if(ndims(image) > 2)
    image_gray = double(rgb2gray(image))/255;
    %image_gray = NormalizeVals(image_gray);
else
    image_gray = double(image)/255;
    %image_gray = NormalizeVals(image_gray);
end

%image_gray = NormalizeVals(image_gray);
figure('Name','Preprocessing: Grayscale');
imshow([image_gray,RefBw]);
title('Preprocessing: Grayscale');

%Displaying the difference between the two
% BwDiff = PercentError(RefBw,image_gray);
% figure('Name','% Difference in Grayscale');
% imagesc(BwDiff*100);
% colorbar;
% title('% Difference in Grayscale');
% BwTotalErr = sum(sum(BwDiff))*100;

%Stage 1: Gaussian Blurring 
image_blurred = imgaussfilt(image_gray, 0.8);

%Displaying the results of blurring
figure('Name','Stage 1:Gaussian Blurring');
imshow([image_blurred,RefBlur]);
title('Stage 1:Gaussian Blurring');

%Displaying the difference between the two
% BlurDiff = PercentError(RefBlur,image_blurred);
% figure('Name','% Difference in Blurring');
% imagesc(BlurDiff*100);
% colorbar;
% title('% Difference in Blurring');
% BlurTotalErr = sum(sum(BlurDiff))*100;


%Stage 2: Calculating Gradients (Without toolbox)
G = fspecial('gaussian',3,0.5); %Generate Gausian Filter
[dx, dy] = gradient(G);         %Calc Gradient of Gausian
Ix = conv2(RefBlur,dx,'same');
Iy = conv2(RefBlur,dy,'same');

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
imagesc([gm,RefMag]);
colorbar;
title('Stage 2a: Gradient Magnitue');
figure('Name','Stage 2b: Gradient Direction');
imagesc([gd,RefTheta]);
colorbar;
title('Stage 2b: Gradient Direction');
tStep2 = toc(tStart);

% MagDiff = PercentError(RefMag,gm);
% figure('Name','% Difference in Magnitude');
% imagesc(MagDiff*100);
% colorbar;
% title('% Difference in Magnitude');
% MagTotalErr = sum(sum(MagDiff))*100;
% 
% ThetaDiff = PercentError(RefTheta,gd);
% figure('Name','% Difference in Theta');
% imagesc(ThetaDiff*100);
% colorbar;
% title('% Difference in Theta');
% BlurTotalErr = sum(sum(ThetaDiff))*100;

%Stage 3 + 4: Edge Extraction / Clustering
image_clusters = CalcEdges(ArraytoList(RefMag),ArraytoList(RefTheta),0.004);
tStep3_4 = toc(tStart) - tStep2

Cluster_Num = unique(image_clusters(:,4)); %Gets each unique cluster
current_num = 1; %holds the offset of the where we're grabbing clusters

    figure('Name','Grouped Segments');
    imshow(image_gray);
    title('Grouped Segments');
    hold on;
    for i = 1:size(Cluster_Num)
        num_of_pts = size(find(image_clusters(:,4) == Cluster_Num(i)),1);
        temp = image_clusters(current_num:num_of_pts+current_num - 1,:);
        plot(temp(:,1),temp(:,2),'*');
        current_num = current_num + num_of_pts; %Add to the offset
    end


%Stage 5: Segmentation 
MinCluster = 4;
FoundSegs   = Segmenter(image_clusters,ArraytoList(RefTheta),ArraytoList(RefMag),image_gray);
tStep5 = toc(tStart) - tStep3_4

%Stage 6: Chain Segments
linked_segments = LinkSegs(FoundSegs);
tStep6 = toc(tStart) - tStep5

%Stage 7: Find Quads
quads = QuadDetection(linked_segments,FoundSegs);
tStep7 = toc(tStart) - tStep6

    %Debug visualization
    figure('Name','Detected Quads with intersections');
    imshow(image_gray);
    title('Detected Quads with intersections');
    hold on;
    for i = 1:size(quads,1)
        Seg1 = [quads(i,1),quads(i,3); quads(i,2), quads(i,4)];
        Seg2 = [quads(i,3),quads(i,5); quads(i,4), quads(i,6)];
        Seg3 = [quads(i,5),quads(i,7); quads(i,6), quads(i,8)];
        Seg4 = [quads(i,7),quads(i,1); quads(i,8), quads(i,2)];
        
        plot(Seg1(1,:),Seg1(2,:),'r-');
        plot(Seg2(1,:),Seg2(2,:),'r-');
        plot(Seg3(1,:),Seg3(2,:),'r-');
        plot(Seg4(1,:),Seg4(2,:),'r-');
        scatter([quads(i,1),quads(i,3),quads(i,5),quads(i,7)],[quads(i,2),quads(i,4),quads(i,6),quads(i,8)],15,'go');
    end

%Stage 8: Decode Quads
Detections = DecodeQuad(quads,RefBw);
tStep8 = toc(tStart) - tStep7

%Stage 9: Remove Duplicates (Skipping For Now)

%Stage 10?: Decode Pose From Detections
Pose = struct('dist',0,'x',0,'y',0,'z',0,'pitch',0,'roll',0,'yaw',0);
Pose = PoseDecoding(Detections,TagSize,Fx,Fy,Px,Py);

tElapsed = toc(tStart)

sprintf('I found %i tag(s)\n',size(Detections));
for NumDet = 1:size(Detections)
    sprintf('Id:%i (Hamming: %i)',Detections(NumDet).id,Detections(NumDet).HD);
    sprintf('distance=%5fm, x=%5f, y=%5f, z=%5f, pitch=%5f, roll=%5f, yaw=%5f',...
        Pose(NumDet).dist,Pose(NumDet).x,Pose(NumDet).y,Pose(NumDet).z,...
        Pose(NumDet).pitch,Pose(NumDet).roll,Pose(NumDet).yaw);
end


function Pose = PoseDecoding(Detections,TagSize,Fx,Fy,Px,Py)
for i = 1:size(Detections)
    TD_getRelativeTandR(TagSize,Fx,Fy,Px,Py,Detections(i))
end
end

function TagPose = TD_getRelativeTandR(TagSize,Fx,Fy,Px,Py,TD_struct)
H = TD_struct.homography;

R20 =  H(3,1);
R21 =  H(3,2);
TZ  =  H(3,3); %Really TZ
R00 = (H(1,1) - Px*R20) / Fx;
R01 = (H(1,2) - Px*R21) / Fx;
TX  = (H(1,3) - Px*TZ)  / Fx;
R10 = (H(2,1) - Py*R20) / Fy;
R11 = (H(2,2) - Py*R21) / Fy;
TY  = (H(2,3) - Py*TZ)  / Fy;

length1 = sqrt(R00^2 + R10^2 + R20^2);
length2 = sqrt(R01^2 + R11^2 + R21^2);
s = (sqrt(length1 * length2))^-1;
if(TZ > 0)
    s = s * -1;
end

R20 = R20 * s;
R21 = R21 * s;
TZ  = TZ  * s;
R00 = R00 * s;
R01 = R01 * s;
TX  = TX  * s;
R10 = R10 * s;
R11 = R11 * s;
TY  = TY  * s;


R02 = R10 * R21 - R20*R11;
R12 = R20 * R01 - R00*R21;
R22 = R00 * R11 - R10*R01;

R = [R00,R01,R02;R10,R11,R12;R20,R21,R22];

[U,~,V] = svd(R);

R = U * V';

R00 = R(1,1);
R01 = R(1,2);
R02 = R(1,3);
R10 = R(2,1);
R11 = R(2,2);
R12 = R(2,3);
R20 = R(3,1);
R21 = R(3,2);
R22 = R(3,3);

TagPose = [R00,R01,R02,TX;R10,R11,R12,TY;R20,R21,R22,TZ;0,0,0,1]
end

%These are helper / utility functions

% function GrayImage = cvtColor(InputImage)
% RedConv   = single(InputImage(:,:,1) *  0.299);
% GreenConv = single(InputImage(:,:,2) *  0.587);
% BlueConv  = single(InputImage(:,:,3) *  0.114);
% 
% GrayImage = RedConv + GreenConv + BlueConv;
% GrayImage = GrayImage / 255;
% end

function output = NormalizeVals(input,Max,Min)
    switch nargin
        case 1
            output = (input-min(input(:)))./(max(input(:))-min(input(:)));
        otherwise
            output = (input-Min)./(Max-Min);
    end
end

function longArray = ArraytoList(Array)
Width = size(Array,2);
Height  = size(Array,1);

longArray = zeros(1,Width*Height);
for i = 1:Height
    StartIdx = ((i-1) * Width)+1;
    EndIdx   = (StartIdx + Width)-1;
    longArray(1,StartIdx:EndIdx) = Array(i,:);
end
end

function Error = PercentError(Correct,Experimental)
difference = Experimental - Correct;
Error = abs(difference./Correct);
end
