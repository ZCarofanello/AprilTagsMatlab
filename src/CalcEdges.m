function Edge = CalcEdges(Magnitude, Direction, MagThr, height, width)
    MinMag = MagThr;
    Edge = [];
    
    FirstEntry = true; %Bool to make sure we don't miss the first entry
    for y = 5:height-5
        for x = 5:width-5
            if(Magnitude(y*width+x) > MinMag)
            %Cost1
             if(Magnitude(y*width+(x+1)) > MinMag)
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction(y*width+(x+1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*width+x;
                    IdB  = y*width+(x+1);
                    Point = [x+1,y];

                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost2
            if(Magnitude((y+1)*width+(x)) > MinMag)
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x);  
                    Point = [x,y+1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost3
            if(Magnitude((y+1)*width+(x+1)) > MinMag)
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x+1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x+1); 
                    Point = [x+1,y+1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost4
            if(Magnitude((y+1)*width+(x-1)) > MinMag && x ~= 2)
                E_Cost = EdgeCost(Direction(y*width+x)...
                    ,Direction((y+1)*width+(x-1)));
                
                %If cost is above threshold then add it to list
                if(E_Cost >= 0)
                	Cost = double(E_Cost);
                    IdA  = y*width+x;
                    IdB  = (y+1)*width+(x-1);  
                    Point = [x-1,y+1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            end
        end
    end
    Edge = sortrows(Edge,1); %Not needed but helps the merger a little
    %Display found Edges
%     figure;
%     imshow(Magnitude);
%     hold on;
%     for i = 1:size(Edge,1)
%         plot(Edge(i,5),Edge(i,4),'r*');
%     end
%     hold off;
    
    Edge = MergeEdges(Edge,Magnitude,Direction); %Merges the detected edges
end

function cost = EdgeCost(Theta0,Theta1)
%Calculates the edge cost between two thetas

maxEdgeCost = (30 * pi) / 180;       %Makes the max edge cost 30 degrees

cost = abs(mod2pi(Theta1 - Theta0)); %The difference between the two edges

if(cost > maxEdgeCost)
    cost = int16(-1); %Essentially makes the edge cost infinite
    return;
end

cost = cost / maxEdgeCost; %The percentage of the cost to the max edge cost

cost = floor(cost * 100); %Scale the cost up by 100
end