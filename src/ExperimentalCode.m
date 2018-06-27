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