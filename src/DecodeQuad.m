function TagDetections = DecodeQuad(quads,GrayImg,Debug)
%Constants to export
blackBorder = 1; dimension = 6;
[height,width] = size(GrayImg);
TagDetections = [];

OC = OpticalCenter(width,height);
for i = 1:size(quads,1)
    
    %To keep me sane
    ThisQuad = [quads(i,1),quads(i,2);quads(i,3),quads(i,4);...
                quads(i,5),quads(i,6);quads(i,7),quads(i,8)];
    
    
    %Initalizing the gray models for this quad
    blackModel = GM_Init();
    whiteModel = GM_Init();
    
    %Initalizing the Homography for this quad
    Quad_H33 = H33_Init(OC);
    Quad_H33 = H33_AddCorrespondence(-1,-1,ThisQuad(1,1),ThisQuad(1,2),Quad_H33);
    Quad_H33 = H33_AddCorrespondence( 1,-1,ThisQuad(2,1),ThisQuad(2,2),Quad_H33);
    Quad_H33 = H33_AddCorrespondence( 1, 1,ThisQuad(3,1),ThisQuad(3,2),Quad_H33);
    Quad_H33 = H33_AddCorrespondence(-1, 1,ThisQuad(4,1),ThisQuad(4,2),Quad_H33);

    dd = 2 * blackBorder + dimension; %Find DD

    %Debug figure
    if(Debug)
        figure('Name','Mapping Points');
        imshow(GrayImg);
        title('Local Mapping points');
        hold on;
    end
    
    for iy = -1:dd
        y = (iy + 0.5) / dd; %Generate local y
        for ix = -1:dd
            x = (ix + 0.5) / dd; %Generate local x
            [px,py,Quad_H33]= Quad_interpolate01(x,y,Quad_H33); %find actual x and y
            
            irx = floor((px) + 0.5); %Get actual x value
            iry = floor((py)+ 0.5); %Get actual y value
            
            %check if it's a valid value
            if(irx < 0 || irx >= width || iry < 0 || iry >= height)
                continue;
            end
            
            v = GrayImg(iry,irx); %Get grayscale pixel value
            
            %Debug visualization
            if(Debug)
                scatter(irx,iry,20,'filled','r');
            end
            
            %if within bounds add observation to either white or black
            if (iy == -1 || iy == dd || ix == -1 || ix == dd)
                whiteModel = GM_addObs(x,y,v,whiteModel);
            elseif (iy == 0 || iy == (dd-1) || ix == 0 || ix == (dd-1))
                blackModel = GM_addObs(x,y,v,blackModel);
            end  
        end
    end
    if (Debug)
    hold off; %debug plot
    end
    
    %Debug figure
    if (Debug)
        figure('Name','Decoding Tag Contents');
        imshow(GrayImg);
        title('Decoding Tag Contents');
        hold on;
    end
    
    bad = false;
    tagCode = uint64(0);
    
    for iy = dimension-1:-1:0
        y = (blackBorder + iy + 0.5)/dd; %Generate local y value
        
        for ix = 0:dimension-1
            x = (blackBorder + ix +0.5) /dd; %Generate local x value
            
            %Find actual x y
            [px,py,Quad_H33] = Quad_interpolate01(x,y,Quad_H33);
            
            irx = floor((px ) + 0.5); %Get actual x value
            iry = floor((py)+ 0.5); %Get actual y value
            
            if( irx < 0 || irx >= width || iry < 0 || iry >= height)
                bad = true;
                continue;
            end
            %Get observations from Black and White Models
            [thrBM,blackModel] = GM_interpolate(x,y,blackModel);
            [thrWM,whiteModel] = GM_interpolate(x,y,whiteModel);
            
            threshold = (thrBM + thrWM) * 0.5; %Create threshold for point
            
            
            v = GrayImg(iry,irx); %Get grayscale value
            
            tagCode = bitshift(tagCode,1); %Shift left to get new value
            if( v > threshold) %if obs is above threshold it's a 1
                tagCode = bitor(tagCode,uint64(1)); 
            end
            
            %Debugging visualization
            if(Debug)
                if( v > threshold )
                    scatter(irx,iry,20,'filled','g');
                else
                    scatter(irx,iry,20,'filled','r');
                end
            end
        end
    end
    if(Debug)
    hold off
    end
    if(~bad)
        TagDetection = TF_Decode(tagCode);
        TagDetection.homography = Quad_H33.H;
        %TagDetection.hxy = Quad_H33.cxy;

        %Correcting the rotation of the tag
%         c = cos(TagDetection.Rotation * (pi/2));
%         s = sin(TagDetection.Rotation * (pi/2));
%         R = zeros(3);
%         R(1,1) = c; R(2,2) = c;
%         R(1,2) = -s;
%         R(2,1) = s;
%         R(3,3) = 1;
% 
%         tmp = TagDetection.homography * R;
%         TagDetection.homography = tmp;

        %Should be the bottom left point of the tag
        [bLx,bLy] = Quad_interpolate01(-1,-1,Quad_H33);
        bestRot = -1;
        bestDist = realmax;
                
        %Add those points to the Output struct
        TagDetection.QuadPts = ThisQuad;
        
        %loop through quad points to see which one is the bottom left pt
        for j = 1:4
            dist = Pt2PtDist(bLx,bLy,ThisQuad(j,1),ThisQuad(j,2));
            if(dist  < bestDist)
                bestDist = dist;
                bestRot = j;
            end
        end
        
        %Kludged way of making the points the correct rotation
        switch bestRot
            case 1
                TagDetection.QuadPts(1,:) = ThisQuad(1,:);
                TagDetection.QuadPts(2,:) = ThisQuad(2,:);
                TagDetection.QuadPts(3,:) = ThisQuad(3,:);
                TagDetection.QuadPts(4,:) = ThisQuad(4,:);
            case 2
                TagDetection.QuadPts(1,:) = ThisQuad(2,:);
                TagDetection.QuadPts(2,:) = ThisQuad(3,:);
                TagDetection.QuadPts(3,:) = ThisQuad(4,:);
                TagDetection.QuadPts(4,:) = ThisQuad(1,:);
            case 3
                TagDetection.QuadPts(1,:) = ThisQuad(3,:);
                TagDetection.QuadPts(2,:) = ThisQuad(4,:);
                TagDetection.QuadPts(3,:) = ThisQuad(1,:);
                TagDetection.QuadPts(4,:) = ThisQuad(2,:);
            case 4
                TagDetection.QuadPts(1,:) = ThisQuad(4,:);
                TagDetection.QuadPts(2,:) = ThisQuad(1,:);
                TagDetection.QuadPts(3,:) = ThisQuad(2,:);
                TagDetection.QuadPts(4,:) = ThisQuad(3,:);
        end

        if(TagDetection.good)
            %Get the center of the tag
            [cx,cy] = Quad_interpolate01(0.5,0.5,Quad_H33);
            %Get the center of the scene
            TagDetection.cxy = [cx,cy];
            %Needed for step 9
            TagDetection.obsPerimeter = [];
            TagDetections = [TagDetections;TagDetection];
        end
    end
end

end

function Homography = H33_Init(OpticalCenter)
Homography = struct('cxy',OpticalCenter,'fA',[],'H',[],'Valid',[]);
Homography.fA = zeros(9); %Fill matrix with zeros
Homography.H  = zeros(3); %Fill matrix with zeros
Homography.Valid = false; %Init variable
end

function H33_struct = H33_AddCorrespondence(Worldx, Worldy, Imagex, Imagey, H33_struct)
H33_struct.Valid = false;

%Make points in relation to optical center
Imagex = Imagex - H33_struct.cxy(1);
Imagey = Imagey - H33_struct.cxy(2);

%From here down is a complicated mess of homography calculations

a03 = -Worldx;
a04 = -Worldy;
a05 = -1;
a06 = Worldx*Imagey;
a07 = Worldy*Imagey;
a08 = Imagey;

H33_struct.fA(4,4) = H33_struct.fA(4,4) + a03*a03;
H33_struct.fA(4,5) = H33_struct.fA(4,5) + a03*a04;
H33_struct.fA(4,6) = H33_struct.fA(4,6) + a03*a05;
H33_struct.fA(4,7) = H33_struct.fA(4,7) + a03*a06;
H33_struct.fA(4,8) = H33_struct.fA(4,8) + a03*a07;
H33_struct.fA(4,9) = H33_struct.fA(4,9) + a03*a08;

H33_struct.fA(5,5) = H33_struct.fA(5,5) + a04*a04;
H33_struct.fA(5,6) = H33_struct.fA(5,6) + a04*a05;
H33_struct.fA(5,7) = H33_struct.fA(5,7) + a04*a06;
H33_struct.fA(5,8) = H33_struct.fA(5,8) + a04*a07;
H33_struct.fA(5,9) = H33_struct.fA(5,9) + a04*a08;
                                          
H33_struct.fA(6,6) = H33_struct.fA(6,6) + a05*a05;
H33_struct.fA(6,7) = H33_struct.fA(6,7) + a05*a06;
H33_struct.fA(6,8) = H33_struct.fA(6,8) + a05*a07;
H33_struct.fA(6,9) = H33_struct.fA(6,9) + a05*a08;
                                          
H33_struct.fA(7,7) = H33_struct.fA(7,7) + a06*a06;
H33_struct.fA(7,8) = H33_struct.fA(7,8) + a06*a07;
H33_struct.fA(7,9) = H33_struct.fA(7,9) + a06*a08;
                                          
H33_struct.fA(8,8) = H33_struct.fA(8,8) + a07*a07;
H33_struct.fA(8,9) = H33_struct.fA(8,9) + a07*a08;
                                          
H33_struct.fA(9,9) = H33_struct.fA(9,9) + a08*a08;

a10 = Worldx;
a11 = Worldy;
a12 = 1;
a16 = -Worldx*Imagex;
a17 = -Worldy*Imagex;
a18 = -Imagex;

H33_struct.fA(1,1) = H33_struct.fA(1,1) + a10*a10;
H33_struct.fA(1,2) = H33_struct.fA(1,2) + a10*a11;
H33_struct.fA(1,3) = H33_struct.fA(1,3) + a10*a12;
H33_struct.fA(1,7) = H33_struct.fA(1,7) + a10*a16;
H33_struct.fA(1,8) = H33_struct.fA(1,8) + a10*a17;
H33_struct.fA(1,9) = H33_struct.fA(1,9) + a10*a18;
                      
H33_struct.fA(2,2) = H33_struct.fA(2,2) + a11*a11;
H33_struct.fA(2,3) = H33_struct.fA(2,3) + a11*a12;
H33_struct.fA(2,7) = H33_struct.fA(2,7) + a11*a16;
H33_struct.fA(2,8) = H33_struct.fA(2,8) + a11*a17;
H33_struct.fA(2,9) = H33_struct.fA(2,9) + a11*a18;
                      
H33_struct.fA(3,3) = H33_struct.fA(3,3) + a12*a12;
H33_struct.fA(3,7) = H33_struct.fA(3,7) + a12*a16;
H33_struct.fA(3,8) = H33_struct.fA(3,8) + a12*a17;
H33_struct.fA(3,9) = H33_struct.fA(3,9) + a12*a18;
                      
H33_struct.fA(7,7) = H33_struct.fA(7,7) + a16*a16;
H33_struct.fA(7,8) = H33_struct.fA(7,8) + a16*a17;
H33_struct.fA(7,9) = H33_struct.fA(7,9) + a16*a18;
                      
H33_struct.fA(8,8) = H33_struct.fA(8,8) + a17*a17;
H33_struct.fA(8,9) = H33_struct.fA(8,9) + a17*a18;
                      
H33_struct.fA(9,9) = H33_struct.fA(9,9) + a18*a18;

a20 = -Worldx*Imagey;
a21 = -Worldy*Imagey;
a22 = -Imagey;
a23 = Worldx*Imagex;
a24 = Worldy*Imagex;
a25 = Imagex;

H33_struct.fA(1,1) = H33_struct.fA(1,1) + a20*a20;
H33_struct.fA(1,2) = H33_struct.fA(1,2) + a20*a21;
H33_struct.fA(1,3) = H33_struct.fA(1,3) + a20*a22;
H33_struct.fA(1,4) = H33_struct.fA(1,4) + a20*a23;
H33_struct.fA(1,5) = H33_struct.fA(1,5) + a20*a24;
H33_struct.fA(1,6) = H33_struct.fA(1,6) + a20*a25;
                                          
H33_struct.fA(2,2) = H33_struct.fA(2,2) + a21*a21;
H33_struct.fA(2,3) = H33_struct.fA(2,3) + a21*a22;
H33_struct.fA(2,4) = H33_struct.fA(2,4) + a21*a23;
H33_struct.fA(2,5) = H33_struct.fA(2,5) + a21*a24;
H33_struct.fA(2,6) = H33_struct.fA(2,6) + a21*a25;
                                          
H33_struct.fA(3,3) = H33_struct.fA(3,3) + a22*a22;
H33_struct.fA(3,4) = H33_struct.fA(3,4) + a22*a23;
H33_struct.fA(3,5) = H33_struct.fA(3,5) + a22*a24;
H33_struct.fA(3,6) = H33_struct.fA(3,6) + a22*a25;
                                          
H33_struct.fA(4,4) = H33_struct.fA(4,4) + a23*a23;
H33_struct.fA(4,5) = H33_struct.fA(4,5) + a23*a24;
H33_struct.fA(4,6) = H33_struct.fA(4,6) + a23*a25;
                                          
H33_struct.fA(5,5) = H33_struct.fA(5,5) + a24*a24;
H33_struct.fA(5,6) = H33_struct.fA(5,6) + a24*a25;
                                          
H33_struct.fA(6,6) = H33_struct.fA(6,6) + a25*a25;

H33_struct.Valid = false;
end

function H33_struct = H33_Compute(H33_struct)
if(H33_struct.Valid)
    return;
end

%Make Matrix Symetric 
[n,m]=size(H33_struct.fA);
B = H33_struct.fA'+H33_struct.fA;
B(1:n+1:end)=diag(H33_struct.fA);
H33_struct.fA = B;

[~,~,V] = svd(H33_struct.fA); %Calc svd of fA

tmp = V(:,size(V,2)); %Get values that we're interested in

%this is to format the values correctly for the next calculations
for i = 1:3
    H33_struct.H(i,:) = tmp((i-1)*3+1:i*3);
end

H33_struct.Valid = true;
end

function [x,y, H33_struct] = H33_Project(Worldx, Worldy, H33_struct)
%This translates the local coordinate system to the world coordinate sys

if(H33_struct.Valid == false)
    %Need an updated calcuation of H
    H33_struct = H33_Compute(H33_struct);
end
x = H33_struct.H(1,1)*Worldx + H33_struct.H(1,2)*Worldy + H33_struct.H(1,3);
y = H33_struct.H(2,1)*Worldx + H33_struct.H(2,2)*Worldy + H33_struct.H(2,3);

z = H33_struct.H(3,1)*Worldx + H33_struct.H(3,2)*Worldy + H33_struct.H(3,3);

x = (x/z) + H33_struct.cxy(1);
y = (y/z) + H33_struct.cxy(2);
end

function [x,y,H33_Struct] = Quad_interpolate01(x,y, H33_Struct)
[x,y,H33_Struct] = H33_Project(2*x-1, 2*y-1,H33_Struct); %Needed for scale correctly
end

function GrayModel = GM_Init()
GrayModel = struct('A',zeros(4),'v',zeros(4,1),'b',zeros(4,1),'nobs',0,'dirty',true);
end

function GrayModel = GM_addObs(x,y,gray,GrayModel)
xy = x*y;
%Homography matrix math for the gray models
GrayModel.A(1,1) = GrayModel.A(1,1) + x*x;
GrayModel.A(1,2) = GrayModel.A(1,2) + x*y;
GrayModel.A(1,3) = GrayModel.A(1,3) + x*xy;
GrayModel.A(1,4) = GrayModel.A(1,4) + x;
GrayModel.A(2,2) = GrayModel.A(2,2) + y*y;
GrayModel.A(2,3) = GrayModel.A(2,3) + y*xy;
GrayModel.A(2,4) = GrayModel.A(2,4) + y;
GrayModel.A(3,3) = GrayModel.A(3,3) + xy*xy;
GrayModel.A(3,4) = GrayModel.A(3,4) + xy;
GrayModel.A(4,4) = GrayModel.A(4,4) + 1;

GrayModel.b(1,1) = GrayModel.b(1,1) + x*gray;
GrayModel.b(2,1) = GrayModel.b(2,1) + y*gray;
GrayModel.b(3,1) = GrayModel.b(3,1) + xy*gray;
GrayModel.b(4,1) = GrayModel.b(4,1) + gray;

GrayModel.nobs = GrayModel.nobs + 1;
GrayModel.dirty = true;
end

function GrayModel = GM_compute(GrayModel)

GrayModel.dirty = false;

%Homography matrix math for the gray model
if(GrayModel.nobs >= 6)
    
    %Make Matrix Symetric 
    [n,m]=size(GrayModel.A);
    B = GrayModel.A'+GrayModel.A;
    B(1:n+1:end)=diag(GrayModel.A);
    GrayModel.A = B;
    
    GrayModel.v = GrayModel.A^-1 * GrayModel.b; %Multiply A^-1 and b
    
else
    GrayModel.v = zeros(4,1);
    GrayModel.v(4,1) = GrayModel.b(3,1) / nobs;
end
end

function [val, GrayModel] = GM_interpolate(x,y,GrayModel)
if(GrayModel.dirty)
    %Need an updated v matrix
    GrayModel = GM_compute(GrayModel);
end

val = GrayModel.v(1)*x + GrayModel.v(2)*y + GrayModel.v(3)*x*y + GrayModel.v(4);
end

function OC = OpticalCenter(height,width)
OC(2) = round(width/2);
OC(1) = round(height/2);
end

function TagDetection = TF_Decode(rCode)
%Constants that need to be exported
errorRecoveryBits = 1;

%Init local variables
bestId = -1;
bestHamming = intmax;
bestRotation = 0;
bestCode = uint64(0);

%Load Tag Family
TagFamily = TagFam36h11();

%Find all the different rotations of this code
rCodes = uint64(zeros(1,4));
rCodes(1) = rCode;
rCodes(2) = rotate90(rCodes(1),6);
rCodes(3) = rotate90(rCodes(2),6);
rCodes(4) = rotate90(rCodes(3),6);

%Search through all the codes and find which ones match the closest to our
%observation
for id = 1:size(TagFamily.Codes,2)
    for rot = 1:4
        thisHamming = HammingDistance(rCodes(rot),TagFamily.Codes(id));
        if(thisHamming < bestHamming)
            bestHamming = thisHamming;
            bestRotation = rot;
            bestId = id-1;
            bestCode = TagFamily.Codes(id);
        end
    end
end

TagDetection = struct('id',[],'HD',[],'Rotation',[],'good',[]...
    ,'obsCode',[],'code',[],'hxy',[],'cxy',[],'QuadPts',[],'homography',[]);

%Export the decoded tag data
TagDetection.id       = bestId;
TagDetection.HD       = bestHamming;
TagDetection.Rotation = bestRotation -1;
TagDetection.good     = (bestHamming <= errorRecoveryBits);
TagDetection.obsCode  = rCode;
TagDetection.code     = bestCode;

end

function RotatedTag = rotate90(w,d)
%Math to 'rotate' the tag id 90 to check for the correct ID

wr = uint64(0);
oneLongLong = uint64(1);

for r = d-1:-1:0
   for c = 0:d-1
       b = r + d*c;
       wr = bitshift(wr,uint64(1));
       if(bitand(uint64(w),bitshift(oneLongLong,uint64(b))) ~= 0)
           wr = bitor(wr,uint64(1));
       end
   end
end
RotatedTag = wr;
end

function Hd = HammingDistance(a,b)
Hd = PopCountReal(bitxor(a,b));
end

function counts = PopCountReal(w)
%Calculates the hamming weight of w
%I needed to spit the unit64 into two unit32s because the largest float
%integer wasn't large enough

%Needed to split the uint64 into two uint32s
TopBitmask = hex2dec('FFFFFFFF00000000');
BotBitmask = hex2dec('00000000FFFFFFFF');

%Snagging top 32 bits
top = bitshift(bitand(w,TopBitmask),-32);
topCount = count(dec2bin(top),'1'); %count all the '1' bits

%Snagging bottom 32 bits
bot = bitand(w,BotBitmask);
botCount = count(dec2bin(bot),'1'); %count all the '1' bits

%Adding the count of bits
counts = topCount + botCount;
end

function [newx,newy] = TD_interpolate(x,y,TD)
% Interpolation for the Tag Detection struct which treats the homography
% differently and has assumptions about the tag rotation

z = TD.homography(3,1)*x + TD.homography(3,2)*y + TD.homography(3,3);
if(z == 0)
    newx = 0;
    newy = 0;
    return;
end
newx = (TD.homography(1,1)*x + TD.homography(1,2)*y + TD.homography(1,3))/z + TD.hxy(1);
newy = (TD.homography(2,1)*x + TD.homography(2,2)*y + TD.homography(2,3))/z + TD.hxy(2);
end

function distance = Pt2PtDist(P1x,P1y,P2x,P2y)
dx = P1x - P2x; %Change in X
dy = P1y - P2y; %Change in Y
distance = sqrt(dx.^2 + dy.^2); %Find the Euclidean distance 
end