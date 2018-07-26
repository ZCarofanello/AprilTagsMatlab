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

PitchPics = [];
for i = 0:60
    PitchPics = [PitchPics;sprintf('../pics/P/%05i.png',i)];
end

PitchObs = [];
StartTime = tic;
for j = 1:size(PitchPics,1)
   CurrentPic = imread(PitchPics(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,0);
   PitchObs = [PitchObs, Pose];
end
ElapsedTime = toc(StartTime);

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

% figure;
% %Pitch Output
% plotYPR(0,PdataC,PitchObs, 0)
% 
% %Roll Output
% plotYPR(1,RdataC,RollObs, 0)
% 
% %Yaw Disp
% plotYPR(2,YdataC,YawObs, 0)
% 
% figure;
% %Pitch Diff
% plotYPR(0,PdataC,PitchObs, 1)
% 
% %Roll Diff
% plotYPR(1,RdataC,RollObs, 1)
% 
% %Yaw Diff
% plotYPR(2,YdataC,YawObs, 1)

function plotYPR(RowNum,CData,MatData,diff)
switch RowNum
case 0
    PlotTitle = 'Pitch Test:';
case 1
    PlotTitle = 'Roll Test:';
case 2
    PlotTitle = 'Yaw Test:';
end
    
if(diff ~= 1)
    subplot(3,3,1+RowNum);
    axis([-30 30 -40 40]);
    if(RowNum == 1)
        line([-30,30],[-30,30],'Color','green')
    else
        line([-30,30],[0,0],'Color','green')
    end
    hold on;
    plot([-30:30],CData(:,6)*(180/pi),'-r');
    plot([-30:30],[MatData(:).pitch]','-b');
    
    title([PlotTitle,'Pitch']);
    legend('Unity','C++','Matlab','location','southeast');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
    
    subplot(3,3,4+RowNum);
    axis([-30 30 -40 40]);
    if(RowNum == 0)
        line([-30,30],[-30,30],'Color','green')
    else
        line([-30,30],[0,0],'Color','green')
    end
    hold on;
    plot([-30:30],CData(:,7)*(180/pi),'-r');
    plot([-30:30],[MatData(:).roll]','-b');
    
    title([PlotTitle,'Roll']);
    legend('Unity','C++','Matlab','location','southeast');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;

    subplot(3,3,7+RowNum);
    axis([-30 30 -40 40]);
    if(RowNum == 2)
        line([-30,30],[-30,30],'Color','green')
    else
        line([-30,30],[0,0],'Color','green')
    end
    hold on;
    plot([-30:30],CData(:,5)*(180/pi),'-r');
    plot([-30:30],[MatData(:).yaw]','-b');
    title([PlotTitle,'Yaw']);
    legend('Unity','C++','Matlab','location','southeast');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
else
    Path = [-30:30];
    
    subplot(3,3,1+RowNum);
    if(RowNum == 1)
        plot([-30:30],((CData(:,6)*(180/pi)) - Path(:)),'-r');
        hold on;
        plot([-30:30],([MatData(:).pitch]' - Path(:)),'-b');
    else
        plot([-30:30],(CData(:,6)*(180/pi)),'-r');
        hold on;
        plot([-30:30],([MatData(:).pitch]'),'-b');
    end

    %axis([-30 30 -10 10]);
    title([PlotTitle,'Pitch Diff']);
    legend('C++','Matlab');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
    
    subplot(3,3,4+RowNum);
    if(RowNum == 0)
        plot([-30:30],(CData(:,7)*(180/pi) - Path(:)),'-r');
        hold on;
        plot([-30:30],([MatData(:).roll]' - Path(:)),'-b');
    else
        plot([-30:30],(CData(:,7)*(180/pi)),'-r');
        hold on;
        plot([-30:30],([MatData(:).roll]'),'-b');
    end
    %axis([-30 30 -10 10]);
    title([PlotTitle,'Roll Diff']);
    legend('C++','Matlab');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;

    subplot(3,3,7+RowNum);
    if(RowNum == 2)
        plot([-30:30],(CData(:,5)*(180/pi) - Path(:)),'-r');
        hold on;
        plot([-30:30],([MatData(:).yaw]' - Path(:)),'-b');
    else
        plot([-30:30],(CData(:,5)*(180/pi)),'-r');
        hold on;
        plot([-30:30],([MatData(:).yaw]'),'-b');
    end
    %axis([-30 30 -10 10]);
    title([PlotTitle,'Yaw Diff']);
    legend('C++','Matlab');
    xlabel('True Rotation (Degrees)') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
end
end
