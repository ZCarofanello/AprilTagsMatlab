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

TestTags = [];
for i = 0:60
    TestTags = [TestTags;sprintf('../pics/P/%05i.png',i)];
end

% for i = 0:60
%     TestTags = [TestTags;sprintf('../pics/R/%05i.png',i)];
% end

for i = 0:60
    TestTags = [TestTags;sprintf('../pics/Y/%05i.png',i)];
end

DemPoses = [];
for j = 1:size(TestTags,1)
   CurrentPic = imread(TestTags(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,1);
   DemPoses = [DemPoses, Pose];
end
