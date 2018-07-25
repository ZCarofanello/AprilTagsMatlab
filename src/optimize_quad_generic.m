function bestQuad = optimize_quad_generic(Family, Image, Quad, Steps)
bestQuad = Quad;
bestScore = QuadGoodness(Family, Image,Quad);

for StepIdx = 1:length(Steps)
    improved = 1;
    
    Max_Repeat = 1;
    
    for repeat = 1:Max_Repeat
        if(improved == 1)
            continue;
        end
        improved = 0;
        
        for i = 1:4
            CurrentStep = Steps(StepIdx);
            nsteps = 1;
            
            this_best_Score = bestScore;
            
            for sx = -nsteps:nsteps
                for sy = -nsteps:nsteps
                    if(sx == 0 && sy == 0)
                        continue;
                    end
                 thisQuad = bestQuad;
                 thisQuad.p(i,1) = bestQuad.p(i,1) + sx*CurrentStep;
                 thisQuad.p(i,2) = bestQuad.p(i,2) + sy*CurrentStep;
                 
                 thisQuad.H33 = H33_Init(thisQuad);
                 
                 if(det(thisQuad.H33.H) == 0)
                     continue;
                 end
                 
                 thisScore = QuadGoodness(Family, Image,thisQuad);
                 
                 if(this_score > this_best_Score)
                     this_best_quad = thisQuad;
                     this_best_score = this_score;
                 end
                end
            end
            
            if(this_best_score > best_score)
                best_quad = this_best_quad;
                best_score = this_best_score;
                improved = 1;
            end
        end
    end
end
                 

end


function Goodness = QuadGoodness(AT_Family, Image, Quad)
white_border = 1;
 dd = 2 * AT_Family.black_border + AT_Family.d; %Find DD

xmin = flintmax; xmax = 0;
ymin = flintmax; ymax = 0;

% figure('Name','Mapping Points');
% imshow(Image);
% title('Local Mapping points');
% hold on;

for iy = -1:dd
    y = (iy + 0.5) / dd; %Generate local y
    for ix = -1:dd
        x = (ix + 0.5) / dd; %Generate local x
        [px,py,Quad.H33]= Quad_interpolate01(x,y,Quad.H33); %find actual x and y

        irx = floor((px) + 0.5); %Get actual x value
        iry = floor((py)+ 0.5); %Get actual y value

        %Debug visualization
        if(0)
            scatter(irx,iry,20,'filled','r');
        end
        
        xmin = min(xmin, irx);
        xmax = max(xmax, irx);
        ymin = min(ymin, iry);
        ymax = max(ymax, iry);
    end
end
if (0)
hold off; %debug plot
end

%Make sure that the points are within bounds
xmin = max(1, xmin);
xmax = min(size(Image,2),xmax);
ymin = max(1, ymin);
ymax = min(size(Image,1),ymax);

W1 = 0; B1 = 0; Wn = 0; Bn = 0;

wsz = dd*white_border;
bsz = dd*AT_Family.black_border;
Hinv = inv(Quad.H33.H);

for y = ymin:ymax
    Hx = Hinv(1,1) * (0.5 + xmin) + Hinv(1,2) + (y + 0.5) + Hinv(1,3);
    Hy = Hinv(2,1) * (0.5 + xmin) + Hinv(2,2) + (y + 0.5) + Hinv(2,3);
    Hh = Hinv(3,1) * (0.5 + xmin) + Hinv(3,2) + (y + 0.5) + Hinv(3,3);

    for x = xmin:xmax
        %project the pixel center.
        tx = 0; ty = 0;
        
        %divide by homogeneous coordinate
        tx = Hx/Hh;
        ty = Hy/Hh;
        
        %if we move x one pixel to the right, here's what
        %happens to our three pre-normalized coordinates.
        Hx = Hx + Hinv(1,1);
        Hy = Hy + Hinv(2,1);
        Hh = Hh + Hinv(3,1);
        
        txa = abs(tx); tya = abs(ty);
        xymax = max(txa, tya);
        
        if(xymax >= 1 + wsz)
            continue;
        end
        
        v = Image(y,x);
        
        %it's within the white border?
        if(xymax >= 1)
            W1 = W1 + v;
            Wn = Wn + 1;
            continue;
        end
        
        %it's within the black border?
        if(xymax >= 1 - bsz)
            B1 = B1 + v;
            Bn = Bn + 1;
            continue;
        end
        
        %Must be a data bit
        continue;
    end
end
%Average Margin between white and black pixels near border
Goodness = 1 * W1 / Wn - 1 * B1 / Bn;
end

function Homography = H33_Init(ThisQuad)
Homography = struct('cxy',[],'fA',[],'H',[],'Valid',[]);

%Initalizing the Homography for this quad
cx = 0; cy = 0;
for k = 1:4
    cx = cx + ThisQuad.p(k,1);
    cy = cy + ThisQuad.p(k,2);
end
cx = cx /4;
cy = cy /4;

Homography.cxy = [cx,cy];
Homography.fA = zeros(9); %Fill matrix with zeros
Homography.H  = zeros(3); %Fill matrix with zeros
Homography.Valid = false; %Init variable

% Homography = H33_AddCorrespondence(-1,-1,ThisQuad.p(1,1),ThisQuad.p(1,2),Homography);
% Homography = H33_AddCorrespondence( 1,-1,ThisQuad.p(2,1),ThisQuad.p(2,2),Homography);
% Homography = H33_AddCorrespondence( 1, 1,ThisQuad.p(3,1),ThisQuad.p(3,2),Homography);
% Homography = H33_AddCorrespondence(-1, 1,ThisQuad.p(4,1),ThisQuad.p(4,2),Homography);

Seg(1) = struct('pt1',ThisQuad.p(1,:),'pt2',ThisQuad.p(2,:),'pts', []);
Seg(2) = struct('pt1',ThisQuad.p(2,:),'pt2',ThisQuad.p(3,:),'pts', []);
Seg(3) = struct('pt1',ThisQuad.p(3,:),'pt2',ThisQuad.p(4,:),'pts', []);
Seg(4) = struct('pt1',ThisQuad.p(4,:),'pt2',ThisQuad.p(1,:),'pts', []);

numOfPts = 256;
Seg(1) = ExtrapolatePts(Seg(1),numOfPts);
Seg(2) = ExtrapolatePts(Seg(2),numOfPts);
Seg(3) = ExtrapolatePts(Seg(3),numOfPts);
Seg(4) = ExtrapolatePts(Seg(4),numOfPts);

PtsInc = 2/numOfPts;
for p = 1:numOfPts
    CurrentInc = p*(PtsInc);
    %Segment 1
    Homography = H33_AddCorrespondence(-1+CurrentInc, -1,Seg(1).pts(p,1),Seg(1).pts(p,2),Homography);
    %Segment 2
    Homography = H33_AddCorrespondence(1, -1+CurrentInc,Seg(2).pts(p,1),Seg(2).pts(p,2),Homography);
    %Segment 3
    Homography = H33_AddCorrespondence(1-CurrentInc,1,Seg(3).pts(p,1),Seg(3).pts(p,2),Homography);
    %Segment 4
    Homography = H33_AddCorrespondence(-1,1-CurrentInc,Seg(4).pts(p,1),Seg(4).pts(p,2),Homography);
end
Homography = H33_AddCorrespondence(0,0,cx,cy,Homography);


end

function LineSeg = ExtrapolatePts(LineSeg, NumOfPts)
Pts = [];
deltax = LineSeg.pt2(1) - LineSeg.pt1(1);
deltay = LineSeg.pt2(2) - LineSeg.pt1(2);

if(deltay == 0)
   inc = deltax / NumOfPts;
   Pts(:,1) = [LineSeg.pt1(1):inc:LineSeg.pt2(1)]';
   Pts(:,2) = LineSeg.pt1(2);
elseif(deltax == 0)
    inc = deltay / NumOfPts;
    Pts(1:NumOfPts,1) = LineSeg.pt1(1);
    Pts(:,2) = [LineSeg.pt1(2):inc:LineSeg.pt2(2)]';
else
    incx = deltax / NumOfPts;
    incy = deltay / NumOfPts;
    Pts(:,1) = [LineSeg.pt1(1):incx:LineSeg.pt2(1)]';
    Pts(:,2) = [LineSeg.pt1(2):incy:LineSeg.pt2(2)]';
end
LineSeg.pts = Pts;
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

H33_struct.fA(4,4) = H33_struct.fA(4,4) + (a03.*a03);
H33_struct.fA(4,5) = H33_struct.fA(4,5) + (a03.*a04);
H33_struct.fA(4,6) = H33_struct.fA(4,6) + (a03.*a05);
H33_struct.fA(4,7) = H33_struct.fA(4,7) + (a03.*a06);
H33_struct.fA(4,8) = H33_struct.fA(4,8) + (a03.*a07);
H33_struct.fA(4,9) = H33_struct.fA(4,9) + (a03.*a08);

H33_struct.fA(5,5) = H33_struct.fA(5,5) + (a04.*a04);
H33_struct.fA(5,6) = H33_struct.fA(5,6) + (a04.*a05);
H33_struct.fA(5,7) = H33_struct.fA(5,7) + (a04.*a06);
H33_struct.fA(5,8) = H33_struct.fA(5,8) + (a04.*a07);
H33_struct.fA(5,9) = H33_struct.fA(5,9) + (a04.*a08);

H33_struct.fA(6,6) = H33_struct.fA(6,6) + (a05.*a05);
H33_struct.fA(6,7) = H33_struct.fA(6,7) + (a05.*a06);
H33_struct.fA(6,8) = H33_struct.fA(6,8) + (a05.*a07);
H33_struct.fA(6,9) = H33_struct.fA(6,9) + (a05.*a08);

H33_struct.fA(7,7) = H33_struct.fA(7,7) + (a06.*a06);
H33_struct.fA(7,8) = H33_struct.fA(7,8) + (a06.*a07);
H33_struct.fA(7,9) = H33_struct.fA(7,9) + (a06.*a08);

H33_struct.fA(8,8) = H33_struct.fA(8,8) + (a07.*a07);
H33_struct.fA(8,9) = H33_struct.fA(8,9) + (a07.*a08);

H33_struct.fA(9,9) = H33_struct.fA(9,9) + (a08.*a08);

a10 = Worldx;
a11 = Worldy;
a12 = 1;
a16 = -Worldx*Imagex;
a17 = -Worldy*Imagex;
a18 = -Imagex;

H33_struct.fA(1,1) = H33_struct.fA(1,1) + (a10.*a10);
H33_struct.fA(1,2) = H33_struct.fA(1,2) + (a10.*a11);
H33_struct.fA(1,3) = H33_struct.fA(1,3) + (a10.*a12);
H33_struct.fA(1,7) = H33_struct.fA(1,7) + (a10.*a16);
H33_struct.fA(1,8) = H33_struct.fA(1,8) + (a10.*a17);
H33_struct.fA(1,9) = H33_struct.fA(1,9) + (a10.*a18);

H33_struct.fA(2,2) = H33_struct.fA(2,2) + (a11.*a11);
H33_struct.fA(2,3) = H33_struct.fA(2,3) + (a11.*a12);
H33_struct.fA(2,7) = H33_struct.fA(2,7) + (a11.*a16);
H33_struct.fA(2,8) = H33_struct.fA(2,8) + (a11.*a17);
H33_struct.fA(2,9) = H33_struct.fA(2,9) + (a11.*a18);

H33_struct.fA(3,3) = H33_struct.fA(3,3) + (a12.*a12);
H33_struct.fA(3,7) = H33_struct.fA(3,7) + (a12.*a16);
H33_struct.fA(3,8) = H33_struct.fA(3,8) + (a12.*a17);
H33_struct.fA(3,9) = H33_struct.fA(3,9) + (a12.*a18);

H33_struct.fA(7,7) = H33_struct.fA(7,7) + (a16.*a16);
H33_struct.fA(7,8) = H33_struct.fA(7,8) + (a16.*a17);
H33_struct.fA(7,9) = H33_struct.fA(7,9) + (a16.*a18);

H33_struct.fA(8,8) = H33_struct.fA(8,8) + (a17.*a17);
H33_struct.fA(8,9) = H33_struct.fA(8,9) + (a17.*a18);

H33_struct.fA(9,9) = H33_struct.fA(9,9) + (a18.*a18);

a20 = -Worldx*Imagey;
a21 = -Worldy*Imagey;
a22 = -Imagey;
a23 = Worldx*Imagex;
a24 = Worldy*Imagex;
a25 = Imagex;

H33_struct.fA(1,1) = H33_struct.fA(1,1) + (a20.*a20);
H33_struct.fA(1,2) = H33_struct.fA(1,2) + (a20.*a21);
H33_struct.fA(1,3) = H33_struct.fA(1,3) + (a20.*a22);
H33_struct.fA(1,4) = H33_struct.fA(1,4) + (a20.*a23);
H33_struct.fA(1,5) = H33_struct.fA(1,5) + (a20.*a24);
H33_struct.fA(1,6) = H33_struct.fA(1,6) + (a20.*a25);

H33_struct.fA(2,2) = H33_struct.fA(2,2) + (a21.*a21);
H33_struct.fA(2,3) = H33_struct.fA(2,3) + (a21.*a22);
H33_struct.fA(2,4) = H33_struct.fA(2,4) + (a21.*a23);
H33_struct.fA(2,5) = H33_struct.fA(2,5) + (a21.*a24);
H33_struct.fA(2,6) = H33_struct.fA(2,6) + (a21.*a25);

H33_struct.fA(3,3) = H33_struct.fA(3,3) + (a22.*a22);
H33_struct.fA(3,4) = H33_struct.fA(3,4) + (a22.*a23);
H33_struct.fA(3,5) = H33_struct.fA(3,5) + (a22.*a24);
H33_struct.fA(3,6) = H33_struct.fA(3,6) + (a22.*a25);

H33_struct.fA(4,4) = H33_struct.fA(4,4) + (a23.*a23);
H33_struct.fA(4,5) = H33_struct.fA(4,5) + (a23.*a24);
H33_struct.fA(4,6) = H33_struct.fA(4,6) + (a23.*a25);

H33_struct.fA(5,5) = H33_struct.fA(5,5) + (a24.*a24);
H33_struct.fA(5,6) = H33_struct.fA(5,6) + (a24.*a25);

H33_struct.fA(6,6) = H33_struct.fA(6,6) + (a25.*a25);

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

[U,~,~] = svd(H33_struct.fA); %Calc svd of fA

tmp = U(:,size(U,2)); %Get values that we're interested in

%this is to format the values correctly for the next calculations
for i = 1:3
    H33_struct.H(i,:) = tmp((i-1)*3+1:i*3);
end

Tx = eye(3);
Tx(1,3) = 0;
Tx(2,3) = 0;

Ty = eye(3);
Ty(1,3) = H33_struct.cxy(1);
Ty(2,3) = H33_struct.cxy(2);

H33_struct.H = Ty * H33_struct.H * Tx;

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

x = (x/z);
y = (y/z);
end

function [x,y,H33_Struct] = Quad_interpolate01(x,y, H33_Struct)
[x,y,H33_Struct] = H33_Project(2*x-1, 2*y-1,H33_Struct); %Needed for scale correctly
end