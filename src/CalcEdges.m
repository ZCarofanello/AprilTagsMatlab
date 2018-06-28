function Edge = CalcEdges(Magnitude, Direction, MagThr)
    Magnitude(Magnitude <= MagThr) = 0;%Makes sure all edges are above threshold
    
    height = size(Magnitude,1);
    width = size(Magnitude,2);
    Edge = [];
    FirstEntry = true;
    for x = 5:height-5
        for y = 5:width-5
            %Cost1
            if(Magnitude(x+1,y) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = y*height+(x+1);
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
            if(Magnitude(x,y+1) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x,y+1));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x);  
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
            if(Magnitude(x+1,y+1) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y+1));
                if(E_Cost >= 0)
                    Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x+1); 
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
            if(Magnitude(x-1,y+1) ~= 0 && x ~= 2)
                E_Cost = EdgeCost(Direction(x,y),Direction(x-1,y+1));
                if(E_Cost >= 0)
                	Cost = double(E_Cost);
                    IdA  = y*height+x;
                    IdB  = (y+1)*height+(x-1);  
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
    Edge = sortrows(Edge,1);
    %Display found Edges
    figure;
    imshow(Magnitude);
    hold on;
    for i = 1:size(Edge,1)
        plot(Edge(i,5),Edge(i,4),'r*');
    end
    hold off;
    
    Edge = MergeEdges(Edge,Magnitude,Direction);
end

function cost = EdgeCost(Theta0,Theta1)
maxEdgeCost = 30 * pi() / 180;
cost = abs(mod2pi(Theta1 - Theta0));
if(cost > maxEdgeCost)
    cost = int16(-1);
    return;
end
cost = single((cost / maxEdgeCost));
cost = int16(cost * 100);
end