function [Pose, Detections] = AprilTag(image,debug)
Debug_Gradient = 0;

if(nargin < 2)
    debug = 0;
end

% if(debug == 1)
%     figure('Name','Original Image');
%     imshow(image);
%     title('Original Image');
% end

%Constants
TagSize = 0.166;
Fx = 600;
Fy = 600;

%Preprocessing to Grayscale
if(ndims(image) > 2)
    image_gray = single(cvtColor(image));
else
    image_gray = single(image);
end

width = size(image_gray,2);
Px = width/2;
height = size(image_gray,1);
Py = height/2;

if(debug == 1)
figure('Name','Preprocessing: Grayscale');
imshow(image_gray);
title('Preprocessing: Grayscale');
end

%Stage 1: Gaussian Blurring (Without toolbox)
G = fspecial('gaussian',3,0.8); %Generate Gausian Filter
image_blurred = conv2(image_gray,G,'same'); %Convolve across image


%Displaying the results of blurring
if(debug == 1)
figure('Name','Stage 1:Gaussian Blurring');
imshow(image_blurred);
title('Stage 1:Gaussian Blurring');
end

%Stage 2: Calculating Gradients (Without toolbox)
dx = [ 0, 0,0;...
       1, 0,-1;...
       0, 0,0];
dy = [ 0, 1,0;...
       0, 0,0;...
       0,-1,0];
Ix = conv2(image_blurred,dx,'same');  %Convolve across x direction of image
Iy = conv2(image_blurred,dy,'same');  %Convolve across y direction of image

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

gm = single(Ix.^2 + Iy.^2);   %Magnitude
gd = single(atan2(Iy,Ix));    %Direction

if(debug == 1)
figure('Name','Stage 2a: Gradient Magnitue');
imagesc(gm);
colorbar;
title('Stage 2a: Gradient Magnitue');
figure('Name','Stage 2b: Gradient Direction');
imagesc(gd);
colorbar;
title('Stage 2b: Gradient Direction');
end

%Stage 3: Edge Extraction
image_edges = CalcEdges(ArraytoList(gm),ArraytoList(gd)...
    ,0.004, size(image,1), size(image,2));
 
image_clusters = MergeEdges(image_edges,ArraytoList(gm),ArraytoList(gd)); %Merges the detected edges
%image_clusters = MergeEdges_mex(image_edges,ArraytoList(gm),ArraytoList(gd)); %Merges the detected edges

if(debug == 1)
%Debug Code for visualization
Cluster_Num = unique(image_clusters(:,4)); %Gets each unique cluster
current_num = 1; %holds the offset of the where we're grabbing clusters

    figure('Name','Grouped Edges');
    imshow(image_gray);
    title('Grouped Edges');
    hold on;
    for i = 1:size(Cluster_Num)
        num_of_pts = size(find(image_clusters(:,4) == Cluster_Num(i)),1);
        temp = image_clusters(current_num:num_of_pts+current_num - 1,:);
        plot(temp(:,1),temp(:,2),'*','LineWidth',2);
        current_num = current_num + num_of_pts; %Add to the offset
    end
end

%Stage 5: Segmentation 
MinCluster = 4;
FoundSegs   = Segmenter(image_clusters,ArraytoList(gd)...
    ,ArraytoList(gm),width,height);

if(debug == 1)
    figure('Name','Segments');
    imshow(image_gray);
    title('Segments');
    hold on;
    %Debug Code
    for k = 1:length(FoundSegs)
        LineColor = [146/255,abs(FoundSegs(k,5))/(4*pi),1];
        plot([FoundSegs(k,1),FoundSegs(k,3)],...
           [FoundSegs(k,2),FoundSegs(k,4)],...
           'LineWidth',2,'color',LineColor);%plot the segment
    end
    hold off;
end
%Stage 6: Chain Segments
linked_segments = LinkSegs(FoundSegs);

%Stage 7: Find Quads
quads = QuadDetection(linked_segments,FoundSegs);

if(debug == 1)
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
        
        plot(Seg1(1,:),Seg1(2,:),'r-','LineWidth',2);
        plot(Seg2(1,:),Seg2(2,:),'r-','LineWidth',2);
        plot(Seg3(1,:),Seg3(2,:),'r-','LineWidth',2);
        plot(Seg4(1,:),Seg4(2,:),'r-','LineWidth',2);
        scatter([quads(i,1),quads(i,3),quads(i,5),quads(i,7)],...
            [quads(i,2),quads(i,4),quads(i,6),quads(i,8)],15,'go');
        scatter([sum(Seg1(1,:))/2,sum(Seg2(1,:))/2,sum(Seg3(1,:))/2,sum(Seg4(1,:))/2],...
            [sum(Seg1(2,:))/2,sum(Seg2(2,:))/2,sum(Seg3(2,:))/2,sum(Seg4(2,:))/2],15,'go');
    end
end

%Stage 8: Decode Quads
Detections = DecodeQuad(quads,image_gray,0);

%Stage 9: Remove Duplicates (Skipping For Now)
%This part checks if the quad points are on top of eachother and then picks
%the detection with the lower hamming distance or the larger one

%Stage 10?: Decode Pose From Detections
Pose = PoseDecoding(Detections,TagSize,Fx,Fy,Px,Py);

if(debug == 1)
sprintf('I found %i tag(s)\n',size(Detections,1))
for NumDet = 1:size(Detections)
    sprintf('Id:%i (Hamming: %i)',Detections(NumDet).id,Detections(NumDet).HD)
    sprintf('distance=%5fm, x=%5f, y=%5f, z=%5f, pitch=%5f, roll=%5f, yaw=%5f',...
        Pose(NumDet).dist,Pose(NumDet).x,Pose(NumDet).y,Pose(NumDet).z,...
        Pose(NumDet).pitch,Pose(NumDet).roll,Pose(NumDet).yaw)
end
end

if(debug == 1)
    %Debug visualization
    figure('Name','Detected Tags');
    imshow(image);
    title('Detected Tags');
    hold on;
    for i = 1:size(Detections)
        plot(Detections(i).QuadPts(1:2,1),Detections(i).QuadPts(1:2,2),'g-','LineWidth',2);
        plot(Detections(i).QuadPts(2:3,1),Detections(i).QuadPts(2:3,2),'r-','LineWidth',2);
        plot(Detections(i).QuadPts(3:4,1),Detections(i).QuadPts(3:4,2),'m-','LineWidth',2);
        plot(Detections(i).QuadPts([4,1],1),Detections(i).QuadPts([4,1],2),'b-','LineWidth',2);
        scatter(Detections(i).cxy(1),Detections(i).cxy(2),100,'r','LineWidth',2);
        text(Detections(i).cxy(1)+10,Detections(i).cxy(2)+5,sprintf('#%i',Detections(i).id),'color','r');
    end
    hold off;
end

end

%These are helper / utility functions

function GrayImage = cvtColor(InputImage)
RedConv   = single(InputImage(:,:,1) *  0.299);
GreenConv = single(InputImage(:,:,2) *  0.587);
BlueConv  = single(InputImage(:,:,3) *  0.114);

GrayImage = RedConv + GreenConv + BlueConv;
GrayImage = GrayImage / 255;
end

function output = NormalizeVals(input,Max,Min)
    switch nargin
        case 1
            output = (input-min(input(:)))./(max(input(:))-min(input(:)));
        otherwise
            output = (input-Min)./(Max-Min);
    end
end

function longArray = ArraytoList(Array)
%Turns a NxM array into a 1xN*M list 
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