clc;
clear;
close all;

bike  = '../pics/bikesgray.jpg';
tag   = '../pics/test_tag.png';
ref   = '../pics/tag_middle.png';
real  = '../pics/real_life_tag.png';
real2 = '../pics/real_life_tag2.jpg';
real3 = '../pics/real_life_tag3.jpg';
prob  = '../pics/00001.png';

PdataC = csvread('../data/pitchTest.dat');
RdataC = csvread('../data/rollTest.dat');
YdataC = csvread('../data/yawTest.dat');

Path = [-30:30]';
TruePitch = [zeros(61,1),zeros(61,1), Path];


figure;
subplot(1,3,1)
plot(PdataC(:,5)*(180/pi),'r-');
hold on;
title('Pitch C++ Data');
plot(PdataC(:,6)*(180/pi),'g-');
plot(PdataC(:,7)*(180/pi),'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,2)
plot(RdataC(:,5)*(180/pi),'r-');
hold on;
title('Roll C++ Data');
plot(RdataC(:,6)*(180/pi),'g-');
plot(RdataC(:,7)*(180/pi),'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,3)
plot(YdataC(:,5)*(180/pi),'r-');
hold on;
title('Yaw C++ Data');
plot(YdataC(:,6)*(180/pi),'g-');
plot(YdataC(:,7)*(180/pi),'b-');
legend('yaw','pitch','roll')
hold off;


figure;
subplot(1,3,1)
plot(PercentError(PdataC(:,5)*(180/pi),0)*100,'r-');
hold on;
title('% diff:Pitch C++');
plot(PercentError(PdataC(:,6)*(180/pi),0)*100,'g-');
plot(PercentError(PdataC(:,7)*(180/pi),Path)*100,'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,2)
plot(PercentError(RdataC(:,5)*(180/pi),0)*100,'r-');
hold on;
title('% diff:Roll C++');
plot(PercentError(RdataC(:,6)*(180/pi),Path)*100,'g-');
plot(PercentError(RdataC(:,7)*(180/pi),0)*100,'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,3)
plot(PercentError((YdataC(:,5)*(180/pi)),Path)*100,'r-');
hold on;
title('% diff:Yaw C++');
plot(PercentError(YdataC(:,6)*(180/pi),0)*100,'g-');
plot(PercentError(YdataC(:,7)*(180/pi),0)*100,'b-');
legend('yaw','pitch','roll')
hold off;


PitchPics = [];
for i = 0:60
    PitchPics = [PitchPics;sprintf('../pics/P/%05i.png',i)];
end

PitchObs = [];
for j = 1:size(PitchPics,1)
   CurrentPic = imread(PitchPics(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,0);
   PitchObs = [PitchObs, Pose];
end

figure;
subplot(1,3,1)
plot([PitchObs(:).yaw]','r-');
hold on;
title('Pitch Matlab Data');
plot([PitchObs(:).pitch]','g-');
plot([PitchObs(:).roll]','b-');
legend('yaw','pitch','roll')
hold off;

RollPics = [];
for i = 0:60
    RollPics = [RollPics;sprintf('../pics/R/%05i.png',i)];
end

RollObs = [];
for j = 1:size(RollPics,1)
   CurrentPic = imread(RollPics(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,0);
   RollObs = [RollObs, Pose];
end

subplot(1,3,2)
plot([RollObs(:).yaw]','r-');
hold on;
title('Roll Matlab Data');
plot([RollObs(:).pitch]','g-');
plot([RollObs(:).roll]','b-');
legend('yaw','pitch','roll')
hold off;

YawPics = [];
for i = 0:60
    YawPics = [YawPics;sprintf('../pics/Y/%05i.png',i)];
end

YawObs = [];
for j = 1:size(YawPics,1)
   CurrentPic = imread(YawPics(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,0);
   YawObs = [YawObs, Pose];
end

subplot(1,3,3)
plot([YawObs(:).yaw]','r-');
hold on;
title('Yaw Matlab Data');
plot([YawObs(:).pitch]','g-');
plot([YawObs(:).roll]','b-');
legend('yaw','pitch','roll')
hold off;

figure;
subplot(1,3,1)
plot([PitchObs(:).yaw]' - 0,'r-');
hold on;
title('diff between actual:Pitch Mat');
plot([PitchObs(:).pitch]' - 0,'g-');
plot([PitchObs(:).roll]' - Path,'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,2)
plot(PercentError([RollObs(:).yaw]',0)*100,'r-');
hold on;
title('diff between actual:Roll Mat');
plot(PercentError([RollObs(:).pitch]',Path)*100,'g-');
plot(PercentError([RollObs(:).roll]',0)*100,'b-');
legend('yaw','pitch','roll')
hold off;

subplot(1,3,3)
plot(PercentError([YawObs(:).yaw]',Path)*100,'r-');
hold on;
title('diff between actual:Yaw Mat');
plot(PercentError([YawObs(:).pitch]', 0)*100,'g-');
plot(PercentError([YawObs(:).roll]',0)*100,'b-');
legend('yaw','pitch','roll')
hold off;

function Error = PercentError(Experimental,Correct)
difference = Experimental - Correct;
Error = abs(difference./Correct);
end
