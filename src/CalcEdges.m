function Edge = CalcEdges(Magnitude, Direction, MagThr)
    %Magnitude(Magnitude <= MagThr) = 0;%Makes sure all edges are above threshold
    MinMag = MagThr;
    height = size(Magnitude,1);
    width = size(Magnitude,2);
    Edge = [];
    FirstEntry = true;
    for x = 5:height-5
        for y = 5:width-5
            if(Magnitude(x,y) > MinMag)
            if(x == 286 && y == 170)
                AHH = 0;
            end
            %Cost1
            %if(Magnitude(y,x+1) ~= 0)
             if(Magnitude(x+1,y) > MinMag)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = y*height+(x+1);
                    Point = [y,x+1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost2
            %if(Magnitude(y+1,x) ~= 0)
            if(Magnitude(x,y+1) > MinMag)
                E_Cost = EdgeCost(Direction(x,y+1),Direction(x,y+1));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x);  
                    Point = [y+1,x];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost3
            %if(Magnitude(y+1,x+1) ~= 0)
            if(Magnitude(x+1,y+1) > MinMag)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y+1));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x+1); 
                    Point = [y+1,x+1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            %Cost4
            %if(Magnitude(y-1,x+1) ~= 0 && y ~= 2)
            if(Magnitude(x-1,y+1) > MinMag && x ~= 2)
                E_Cost = EdgeCost(Direction(x,y),Direction(x-1,y+1));
                if(E_Cost >= 0)
                	Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x-1);  
                    Point = [y+1,x-1];
                    
                    if(FirstEntry)
                        Edge = [Cost,IdA,IdB,Point];
                        FirstEntry = false;
                    else
                        Edge = [Edge;Cost,IdA,IdB,Point];
                    end
                end
            end
            end
            if(length(Edge) >= 21208)
                AHHH = 0;
            end
        end
    end
    Edge = sortrows(Edge,1);
    %Display found Edges
%     figure;
%     imshow(Magnitude);
%     hold on;
%     for i = 1:size(Edge,1)
%         plot(Edge(i,5),Edge(i,4),'r*');
%     end
%     hold off;
    
    Edge = MergeEdges(Edge,Magnitude,Direction);
end

function cost = EdgeCost(Theta0,Theta1)
maxEdgeCost = (30 * pi) / 180;
cost = abs(mod2pi(Theta1 - Theta0));
if(cost > maxEdgeCost)
    cost = int16(-1);
    return;
end
cost = cost / maxEdgeCost;
cost = floor(cost * 100);
end