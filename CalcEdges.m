function Edge = CalcEdges(Magnitude, Direction)
    height = size(Magnitude,1);
    width = size(Magnitude,2);
    Edge = struct('IdA',[],'IdB',[],'Point',[],'Cost',[]);
    TempStruct = struct('IdA',[],'IdB',[],'Point',[],'Cost',[]);
    FirstEntry = true;
    for x = 2:height-2
        for y = 2:width-2
            %Cost1
            if(Magnitude(x+1,y) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y));
                if(E_Cost >= 0)
                    if(FirstEntry)
                        Edge.Cost = E_Cost;
                        %Edge.IdA  = y+x*height;
                        Edge.IdA  = y*width+x;
                        %Edge.IdB  = y+(x+1)*height;
                        Edge.IdB  = y*width+(x+1); 
                        Edge.Point = [x+1,y];
                        FirstEntry = false;
                    else
                        TempStruct.Cost = E_Cost;
                        %TempStruct.IdA  = y+x*height;
                        TempStruct.IdA  = y*width+x;
                        %TempStruct.IdB  = y+(x+1)*height;
                        TempStruct.IdB  = y*width+(x+1);  
                        TempStruct.Point = [x+1,y];
                        Edge = [Edge,TempStruct];
                    end
                    continue;
                end
            end
            %Cost2
            if(Magnitude(x,y+1) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x,y+1));
                if(E_Cost >= 0)
                    if(FirstEntry)
                        Edge.Cost = E_Cost;
                        %Edge.IdA  = (y)+x*height;
                        Edge.IdA  = y*width+x;
                        %Edge.IdB  = (y+1)+(x)*height;
                        Edge.IdB  = (y+1)*width+(x); 
                        Edge.Point = [x,y+1];
                        FirstEntry = false;
                    else
                        TempStruct.Cost = E_Cost;
                        %TempStruct.IdA  = (y)+x*height;
                        TempStruct.IdA  = y*width+x;
                        %TempStruct.IdB  = (y+1)+(x)*height;
                        TempStruct.IdB  = (y+1)*width+(x);  
                        TempStruct.Point = [x,y+1];
                        Edge = [Edge,TempStruct];
                    end
                    continue;
                end
            end
            %Cost3
            if(Magnitude(x+1,y+1) ~= 0)
                E_Cost = EdgeCost(Direction(x,y),Direction(x+1,y+1));
                if(E_Cost >= 0)
                    if(FirstEntry)
                        Edge.Cost = E_Cost;
                        %Edge.IdA  = (y)+x*height;
                        Edge.IdA  = y*width+x;
                        %Edge.IdB  = (y+1)+(x+1)*height;
                        Edge.IdB  = (y+1)*width+(x+1); 
                        Edge.Point = [x+1,y+1];
                        FirstEntry = false;
                    else
                        TempStruct.Cost = E_Cost;
                        %TempStruct.IdA  = (y)+x*height;
                        TempStruct.IdA  = y*width+x;
                        %TempStruct.IdB  = (y+1)+(x+1)*height;
                        TempStruct.IdB  = (y+1)*width+(x+1); 
                        TempStruct.Point = [x+1,y+1];
                        Edge = [Edge,TempStruct];
                    end
                    continue;
                end
            end
            %Cost4
            if(Magnitude(x-1,y+1) ~= 0 && x ~= 2)
                E_Cost = EdgeCost(Direction(x,y),Direction(x-1,y+1));
                if(E_Cost >= 0)
                    if(FirstEntry)
                        Edge.Cost = E_Cost;
                        %Edge.IdA  = (y)+x*height;
                        Edge.IdA  = y*width+x;
                        %Edge.IdB  = (y+1)+(x-1)*height;
                        Edge.IdB  = (y+1)*width+(x-1);  
                        Edge.Point = [x-1,y+1];
                        FirstEntry = false;
                    else
                        TempStruct.Cost = E_Cost;
                        %TempStruct.IdA  = (y)+x*height;
                        TempStruct.IdA  = y*width+x;
                        %TempStruct.IdB  = (y+1)+(x-1)*height;
                        TempStruct.IdB  = (y+1)*width+(x-1); 
                        TempStruct.Point = [x-1,y+1];
                        Edge = [Edge,TempStruct];
                    end
                    continue;
                end
            end
        end
    end
    %Edge = sortrows(Edge,4);
    %Display found Edges
    figure;
    imshow(Magnitude);
    hold on;
    for i = 1:size(Edge,2)
        plot(Edge(i).Point(2),Edge(i).Point(1),'r*');
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

function value = mod2pi(value, ref)
if nargin == 1
    value = value - 2*pi*floor( (value+pi)/(2*pi) );
    return;
end
shifted = value - ref;
value = (shifted - 2*pi*floor( (shifted+pi)/(2*pi) ))+ref;
end