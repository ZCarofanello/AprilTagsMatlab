%Calculating Cost using matrix math instead of for loops
% function Edge = CalcEdges(Magnitude, Direction)
%     VertBuf = zeros(size(Magnitude,1),1);
%     HorzBuf = zeros(1,size(Magnitude,2));
%     Cost = zeros(size(Magnitude,1),size(Magnitude,2),4);
%     
%     CostMag = horzcat(VertBuf,Magnitude(:,2:size(Magnitude,2))); %Vert
%     CostDir = horzcat(VertBuf,Direction(:,2:size(Magnitude,2))); %Vert
%     Cost(:,:,1) = arrayfun(@EdgeCost,Direction,CostDir,CostMag);
%     Cost1Debug = Cost(:,:,1);
%     figure;
%     image(Cost1Debug);
%     
%     CostMag = vertcat(HorzBuf,Magnitude(2:size(Magnitude,1),:)); %Horz
%     CostDir = vertcat(HorzBuf,Direction(2:size(Magnitude,1),:)); %Horz
%     Cost(:,:,2) = arrayfun(@EdgeCost,Direction,CostDir,CostMag);
%     Cost2Debug = Cost(:,:,2);
%     figure;
%     image(Cost2Debug);
%     
%     CostMag = vertcat(HorzBuf,Magnitude(2:size(Magnitude,1),:)); %Horz
%     CostMag = horzcat(VertBuf,CostMag(:,2:size(Magnitude,2)));   %Vert
%     
%     CostDir = vertcat(HorzBuf,Direction(2:size(Magnitude,1),:)); %Horz
%     CostDir = horzcat(VertBuf,CostDir(:,2:size(Magnitude,2)));   %Vert
%     
%     Cost(:,:,3) = arrayfun(@EdgeCost,Direction,CostDir,CostMag);
%     Cost3Debug = Cost(:,:,3);
%     figure;
%     image(Cost3Debug);
%     
%     CostMag = vertcat(HorzBuf,Magnitude(2:size(Magnitude,1),:)); %Horz
%     CostMag = horzcat(CostMag(:,1:size(Magnitude,2)-1),VertBuf);   %Vert
%     
%     CostDir = vertcat(HorzBuf,Direction(2:size(Magnitude,1),:)); %Horz
%     CostDir = horzcat(CostDir(:,1:size(Magnitude,2)-1),VertBuf);   %Vert
%     
%     Cost(:,:,4) = arrayfun(@EdgeCost,Direction,CostDir,CostMag);
%     Cost4Debug = Cost(:,:,4);
%     figure;
%     image(Cost4Debug);
% 
% end

%Different Angle estimation code 
% if(abs(R20) ~= 1)
%     theta1 = -asin(R20);
%     theta2 = pi - theta1;
%     psi1 = mod2pi(atan2(R21/cos(theta1),R22/cos(theta1)));
%     psi2 = mod2pi(atan2(R21/cos(theta2),R22/cos(theta2)));
%     phi1 = mod2pi(atan2(R10/cos(theta1),R00/cos(theta1)));
%     phi2 = mod2pi(atan2(R10/cos(theta2),R10/cos(theta2)));
% else
%     phi1 = 0;
%     if(R20 == -1)
%         theta1 = pi/2;
%         psi1 = mod2pi(atan2(R01,R02));
%     else
%         theta1 = -pi/2;
%         psi1 = mod2pi(atan2(-R01,-R02));
%     end
%     theta2 = 0; psi2 = 0; phi2 = 0;
% end
% solution1 = [theta1*(180/pi),psi1*(180/pi),phi1*(180/pi)];
% solution2 = [theta2*(180/pi),psi2*(180/pi),phi2*(180/pi)];