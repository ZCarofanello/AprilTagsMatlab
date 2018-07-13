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
for i = 0:50
    TestTags = [TestTags;sprintf('../pics/tags/%05i.png',i)];
end

error = 0;
for j = 1:size(TestTags,1)
   CurrentPic = imread(TestTags(j,:));
   [Pose,Detection] = AprilTag(CurrentPic,0);
   if (Detection.id ~= (j - 1))
       error = error + 1;
   end
end
