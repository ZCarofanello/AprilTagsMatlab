function Edge = CalcEdges(Magnitude, Direction, MagThr, height, width)
    MinMag = MagThr;
    Edge = nan(50000,5);
    EdgeCnt = 1;
    MagLogical = Magnitude > MinMag;
    
    for y = 5:height-5
        for x = 5:width-5
            if(MagLogical(y*width+x))
            %Cost1
             if(MagLogical(y*width+(x+1)))
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction(y*width+(x+1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = E_Cost;
                    IdA  = y*width+x;
                    IdB  = y*width+(x+1);
                    Point = [x+1,y];
                    Edge(EdgeCnt,:) = [Cost,IdA,IdB,Point];
                    EdgeCnt = EdgeCnt + 1;
                end
            end
            %Cost2
            if(MagLogical((y+1)*width+(x)))
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = E_Cost;
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x);  
                    Point = [x,y+1];
                    Edge(EdgeCnt,:) = [Cost,IdA,IdB,Point];
                    EdgeCnt = EdgeCnt + 1;
                end
            end
            %Cost3
            if(MagLogical((y+1)*width+(x+1)))
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x+1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = E_Cost;
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x+1); 
                    Point = [x+1,y+1];
                    Edge(EdgeCnt,:) = [Cost,IdA,IdB,Point];
                    EdgeCnt = EdgeCnt + 1;
                end
            end
            %Cost4
            if(MagLogical((y+1)*width+(x-1)) && x ~= 2)
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x-1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                	Cost = E_Cost;
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x-1);  
                    Point = [x-1,y+1];
                    Edge(EdgeCnt,:) = [Cost,IdA,IdB,Point];
                    EdgeCnt = EdgeCnt + 1;
                end
            end
            end
        end
    end
    Edge(isnan(Edge(:,2)),:) = [];
    
    Edge = sortrows(Edge,1); %Not needed but helps the merger a little
    %Display found Edges
%     figure;
%     imshow(gray_image);
%     hold on;
%     scatter(Edge(:,4),Edge(:,5),'r*');
%     hold off;
    
end

function cost = EdgeCost(Theta0,Theta1)
%Calculates the edge cost between two thetas

maxEdgeCost = (30 * pi) / 180;       %Makes the max edge cost 30 degrees

cost = abs(mod2pi(Theta1 - Theta0)); %The difference between the two edges

if(cost > maxEdgeCost)
    cost = (-1); %Essentially makes the edge cost infinite
    return;
end

cost = cost / maxEdgeCost; %The percentage of the cost to the max edge cost

cost = floor(cost * 100); %Scale the cost up by 100
end