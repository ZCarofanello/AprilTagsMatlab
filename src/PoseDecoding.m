function Pose = PoseDecoding(Detections,TagSize,Fx,Fy,Px,Py)
Pose = [];
for i = 1:size(Detections)
    %Loop through all the valid detections
    Pose = [Pose;TD_getRelativeTandR(TagSize,Fx,Fy,Px,Py,Detections(i))];
end
end

function TagPose = TD_getRelativeTandR(TagSize,Fx,Fy,Px,Py,TD_struct)
Pose = struct('dist',0,'x',0,'y',0,'z',0,'pitch',0,'roll',0,'yaw',0);

CameraMatrix = [Fx, 0, Px;...
                0 ,Fy, Py;...
                0 , 0,  1];

%Fx = -Fx;
H = TD_struct.homography.H;

%Using the homography matrix to construct a rotation and translation matrix
R20 =  H(3,1);
R21 =  H(3,2);
TZ  =  H(3,3);
R00 = (H(1,1) - Px*R20) / Fx;
R01 = (H(1,2) - Px*R21) / Fx;
TX  = (H(1,3) - Px*TZ)  / Fx;
R10 = (H(2,1) - Py*R20) / Fy;
R11 = (H(2,2) - Py*R21) / Fy;
TY  = (H(2,3) - Py*TZ)  / Fy;

%compute the scale by requiring that the rotation columns are unit length
%(Use geometric average of the two length vectors we have)
length1 = sqrt(R00^2 + R10^2 + R20^2);
length2 = sqrt(R01^2 + R11^2 + R21^2);
s = 1/(sqrt(length1 * length2));

%get sign of S by requiring the tag to be in front the camera;
%we assume camera looks in the -Z direction. 
if(TZ > 0)
    s = s * -1;
end

%Normalizing our coefficents by the scaling factor
R20 = R20 * s;
R21 = R21 * s;
TZ  = TZ  * s;
R00 = R00 * s;
R01 = R01 * s;
TX  = TX  * s;
R10 = R10 * s;
R11 = R11 * s;
TY  = TY  * s;

%now recover [R02 R12 R22] by noting that it is the cross product of the other two columns.
R02 = R10 * R21 - R20*R11;
R12 = R20 * R01 - R00*R21;
R22 = R00 * R11 - R10*R01;

%Improve rotation matrix by applying polar decomposition.
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
%End of polar decomposition

%Attempting to account for scale of tag
scale = TagSize/2;
TX  = TX * scale;
TY  = TY * scale;
TZ  = TZ * scale;

%Finished Rotation and translation matrix
RT = [R00,R01,R02,TX;R10,R11,R12,TY;R20,R21,R22,TZ;0,0,0,1];
R = [R00,R01,R02;R10,R11,R12;R20,R21,R22];

RotationFix = [1,0,0;0,-1,0;0,0,1];
FrameFix = [0,0,1;-1,0,0;0,0,-1];
%R = R * FrameFix;
R = R * RotationFix;

[yaw_t,pitch_t,roll_t] = rotation2Euler(R);

%Outputting Calculated Pose
Pose.yaw   = -StandardRad(yaw_t + pi)*(180/pi);
Pose.pitch = -(pitch_t)*(180/pi);
Pose.roll  = StandardRad(roll_t)*(180/pi);
Pose.x     = TX;
Pose.y     = TY;
Pose.z     = TZ;
Pose.dist = sqrt(TX^2 + TY^2 + TZ^2);

TagPose = Pose;
end

function angle = StandardRad(angle)
if(angle >= 0)
    angle = mod(angle+pi, 2*pi) - pi;
else
    angle = mod(angle-pi, -2*pi) + pi;
end
end

function [yaw,pitch,roll] = rotation2Euler(wRo)
yaw = StandardRad(atan2(wRo(2,1),wRo(1,1)));
c = cos(yaw);
s = sin(yaw);
pitch = StandardRad(atan2(-wRo(3,1), wRo(1,1)*c + wRo(2,1)*s));
roll  = StandardRad(atan2(wRo(1,3)*s - wRo(2,3)*c, -wRo(1,2)*s + wRo(2,2)*c));
end