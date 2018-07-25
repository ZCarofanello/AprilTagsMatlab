clear;
close all;
clc;

load('../output/AllData.mat');
PdataC = csvread('../data/pitchTest.dat');
RdataC = csvread('../data/rollTest.dat');
YdataC = csvread('../data/yawTest.dat');


figure;

%Pitch Output
plotYPR(0,PdataC,PitchObs, 0)

%Roll Output
plotYPR(1,RdataC,RollObs, 0)

%Yaw Disp
plotYPR(2,YdataC,YawObs, 0)

figure;
%Pitch Diff
plotYPR(0,PdataC,PitchObs, 1)

%Roll Diff
plotYPR(1,RdataC,RollObs, 1)

%Yaw Diff
plotYPR(2,YdataC,YawObs, 1)


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
    axis([-30 30 -30 30]);
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
    axis([-30 30 -30 30]);
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
    axis([-30 30 -30 30]);
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

function plotXYZ(RowNum,CData,MatData,diff)
switch RowNum
case 0
    PlotTitle = 'X Test:';
case 1
    PlotTitle = 'Y Test:';
case 2
    PlotTitle = 'Z Test:';
end
    
if(diff ~= 1)
    subplot(3,3,1+RowNum);
    if(RowNum == 1)
        line([1,61],[-30,30],'Color','green')
    else
        line([1,61],[0,0],'Color','green')
    end
    hold on;
    plot(CData(:,6)*(180/pi),'-r');
    plot([MatData(:).pitch]','-b');
    axis([1 61 -30 30]);
    title([PlotTitle,'Pitch']);
    legend('Unity','C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
    
    subplot(3,3,4+RowNum);
    if(RowNum == 0)
        line([1,61],[-30,30],'Color','green')
    else
        line([1,61],[0,0],'Color','green')
    end
    hold on;
    plot(CData(:,7)*(180/pi),'-r');
    plot([MatData(:).roll]','-b');
    axis([1 61 -40 40]);
    title([PlotTitle,'Roll']);
    legend('Unity','C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;

    subplot(3,3,7+RowNum);
    if(RowNum == 2)
        line([1,61],[-30,30],'Color','green')
    else
        line([1,61],[0,0],'Color','green')
    end
    hold on;
    plot(CData(:,5)*(180/pi),'-r');
    plot([MatData(:).yaw]','-b');
    axis([1 61 -40 40]);
    title([PlotTitle,'Yaw']);
    legend('Unity','C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
else
    Path = [-30:30];
    
    subplot(3,3,1+RowNum);
    if(RowNum == 1)
        plot(((CData(:,6)*(180/pi)) - Path(:)),'-r');
        hold on;
        plot(([MatData(:).pitch]' - Path(:)),'-b');
    else
        plot((CData(:,6)*(180/pi)),'-r');
        hold on;
        plot(([MatData(:).pitch]'),'-b');
    end

    axis([1 61 -10 10]);
    title([PlotTitle,'Pitch Diff']);
    legend('C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
    
    subplot(3,3,4+RowNum);
    if(RowNum == 0)
        plot((CData(:,7)*(180/pi) - Path(:)),'-r');
        hold on;
        plot(([MatData(:).roll]' - Path(:)),'-b');
    else
        plot((CData(:,7)*(180/pi)),'-r');
        hold on;
        plot(([MatData(:).roll]'),'-b');
    end
    axis([1 61 -10 10]);
    title([PlotTitle,'Roll Diff']);
    legend('C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;

    subplot(3,3,7+RowNum);
    if(RowNum == 2)
        plot((CData(:,5)*(180/pi) - Path(:)),'-r');
        hold on;
        plot(([MatData(:).yaw]' - Path(:)),'-b');
    else
        plot((CData(:,5)*(180/pi)),'-r');
        hold on;
        plot(([MatData(:).yaw]'),'-b');
    end
    axis([1 61 -10 10]);
    title([PlotTitle,'Yaw Diff']);
    legend('C++','Matlab');
    xlabel('Picture #') % x-axis label
    ylabel('degrees') % y-axis label
    hold off;
end
end

