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
    
    SegLength = Pt2PtDist(LineTemp(1,1),LineTemp(1,2),LineTemp(1,3),LineTemp(1,4));
    
    if(MinDist < SegLength)
        LineTemp(6) = SegLength;
        LineTemp = FindDirection(LineTemp,temp,Theta);
        segments = [segments;LineTemp]; %Add to the good segments
        LineColor = abs(LineTemp(5))/(2*pi);
        plot([LineTemp(1),LineTemp(3)],[LineTemp(2),LineTemp(4)],'Color',[146/255,LineColor,1]); %plot the segment
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
      %|  X1  |  Y1  |  X2  |  Y2  | theta | length |
line = [ min_x, min_y, max_x, max_y,0,0];
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
dx = P1x - P2x; %Change in X
dy = P1y - P2y; %Change in Y
distance = sqrt(dx^2 + dy^2); %Find the Euclidean distance 
end

function LineTemp = FindDirection(LineTemp,temp,gd)
dx = LineTemp(3) - LineTemp(1); %Find the change in x
dy = LineTemp(4) - LineTemp(2); %Find the change in y

tmpTheta = atan2(dy,dx);        %Temp direction of the line
theta = gd(temp(:,1),temp(:,2));%Get all the thetas of the line

err = mod2pi(theta - tmpTheta); %Calculate the error of our assumed direction

%Variables for our votes
noflip = 0;
flip = 0;
for i = 1:size(temp)
   if(err(i) < 0) %If the error is negative vote for no flip
       noflip = noflip + temp(i,3);
   else           %If the error is positive vote for to flip
       flip = flip + temp(i,3);
   end
end

if (flip > noflip) %If it's flipped add PI
    LineTemp(5) = tmpTheta + pi;
else
    LineTemp(5) = tmpTheta;
end

dot = dx*cos(LineTemp(5)) + dy*sin(LineTemp(5)); %Check if it's the right dir

if(dot > 0) %If not flip it
    tmpX = LineTemp(1);
    LineTemp(1) = LineTemp(3);
    LineTemp(3) = tmpX;
    
    tmpY = LineTemp(2);
    LineTemp(2) = LineTemp(4);
    LineTemp(4) = tmpY;
end

end