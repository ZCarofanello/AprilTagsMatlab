function Clusters = MergeEdges(Edges,Magnitude,Direction)
    %Constants to export sometime
    thetaThr = 50;
    magThr = 1200;
    
    %Get the width and height of the iamge
    width = size(Magnitude,2);
    height = size(Magnitude,1);
    
    %Reshape the Magnitude and Directions of the arrays
    tmin = reshape(Direction,1,[]);
    tmax = reshape(Direction,1,[]);
    mmin = reshape(Magnitude,1,[]);
    mmax = reshape(Magnitude,1,[]);
    
    %Create the unionfind vector which is pre allocated for speed
    SimpleUF = [(1:width*height); ones(1,width*height);(1:width*height)];
    
    for i = 1:size(Edges,2)
        ida = Edges(i).IdA;
        idb = Edges(i).IdB;
        
        [ida, SimpleUF] = getRepresentative(SimpleUF,ida); %gets rep
        [idb, SimpleUF] = getRepresentative(SimpleUF,idb); %gets rep
        
        if(ida == idb) %It's already connected!
            continue;
        end
        
        sza = SimpleUF(2,ida); %Get the size of tree a
        szb = SimpleUF(2,idb); %Get the size of tree b
        
        tmina = tmin(ida); tmaxa = tmax(ida); %finds the max and mins of a
        tminb = tmin(idb); tmaxb = tmax(idb); %finds the max and mins of b
        
        costa = tmaxa-tmina; %Intermediate cost value
        costb = tmaxb-tminb; %Intermediate cost value
        
        bshift = mod2pi((tmina+tmaxa)/2, (tminb+tmaxb)/2) - (tminb+tmaxb)/2; %Makes sure that the angles aren't more than +/- PI
        
        tminab = min(tmina, tminb+bshift); %Theta min
        tmaxab = max(tmaxa, tmaxb+bshift); %Theta max
        
        mminab = min(mmin(ida), mmin(idb)); %Mag min
        mmaxab = max(mmax(ida), mmax(idb)); %Mag max
        
        costab = (tmaxab - tminab); %Intermediate cost value for a and b
        
        %Magic Values that I need to understand more :)
        Value1 = (costab <= min(costa,costb) + thetaThr/(sza+szb));
        Value2 = (mmaxab-mminab) <= min(mmax(ida)-mmin(ida), mmax(idb)-mmin(idb)) + magThr/(sza+szb);
        
        if(Value1 && Value2)
            [SimpleUF, idab] = connectNodes(SimpleUF, ida, idb);
            
            tmin(idab) = tminab;
            tmax(idab) = tmaxab;
            
            mmin(idab) = mminab;
            mmin(idab) = mmaxab;
        end
    end
    %Not terribly accurate because this is taking simple difference
    %VisualizeClusters(SimpleUF,width,height);
    
    %Export the clusters
    Clusters = ExportClusters(SimpleUF,Magnitude, Edges);
end

% Gets the representative of the node
function [root,UpdatedArray] = getRepresentative(UFArray,NodeId)
    if(UFArray(1,NodeId) == NodeId) %If it is it's own rep return
        root = NodeId;          %No changes
        UpdatedArray = UFArray; %Send back the array
        return;
    end
root = getRepresentative(UFArray,UFArray(1,NodeId)); %Recurse
UFArray(1,NodeId) = root; %Flatten the tree
UpdatedArray = UFArray;   %Return the updated array
end

%connects and merges the two trees together
function [NewUFArray,root] = connectNodes(UFArray, aId,bId)

    [aRoot, UFArray] = getRepresentative(UFArray,aId); %Get rep of a
    [bRoot, UFArray] = getRepresentative(UFArray,bId); %Get rep of b
    NewUFArray = UFArray; %copy the array

    if(aRoot==bRoot) %It's already connected!
        root=aRoot;  %Return the root
        return;
    end
    
    if(UFArray(2,aRoot) > UFArray(2,bRoot)) %Larger tree wins!
        NewUFArray(1,bRoot) = aRoot; %Set the new root
        NewUFArray(2,aRoot) = NewUFArray(2,aRoot) + NewUFArray(2,bRoot); %Add the sizes together
        root=aRoot; %Return the new root
        return;
     else
        NewUFArray(1,aRoot) = bRoot; %Set the new root
        NewUFArray(2,bRoot) = NewUFArray(2,aRoot) + NewUFArray(2,bRoot); %Add the sizes together
        root=bRoot; %Return the new root
        return;
    end
end

function value = mod2pi(value, ref)
if nargin == 1
    value = value - 2*pi*floor( (value+pi)/(2*pi) );
    return;
end
shifted = value - ref;
value = (shifted - 2*pi*floor( (shifted+pi)/(2*pi) ))+ref;
end

%Formatting the clusters as a list with points to make it easier later
function ClusterList = ExportClusters(UF_Array,Magnitude,Edges)
    %Need to export these constants
    MinCluster = 4;

    %find clusters that have more than the MinSeg
    Valid_Clusters = UF_Array(1,(UF_Array(2,:) >= MinCluster));

    %Create a logical array for faster indexing / display
    logical_arr = ismember(UF_Array(1,:),Valid_Clusters);
    
    ClusterList = [];
    FirstEntry = true;
    for i = 1:size(Edges,2)-1
        if(logical_arr(Edges(i).IdA))
            EdgeCluster = UF_Array(1,Edges(i).IdA);
            EdgeMag = Magnitude(Edges(i).Point(1),Edges(i).Point(2));
            EdgeX = Edges(i).Point(2);
            EdgeY = Edges(i).Point(1);
            if(FirstEntry)
                ClusterList = [EdgeX,EdgeY,EdgeMag,EdgeCluster];
                FirstEntry = false;
            else
                ClusterList = [ClusterList;EdgeX,EdgeY,EdgeMag,EdgeCluster];
            end
        end
    end
    figure;
    imshow(Magnitude);
    hold on;
    plot(ClusterList(:,1),ClusterList(:,2),'r*');
    hold off;
    ClusterList = sortrows(ClusterList,4);
end

% %This is a messy function to verify that it's spitting out valid data
% function VisualizeClusters(UF_Array, Width,Height)
%     VertBuf = zeros(Height,1);
%     HorzBuf = zeros(1,Width);
%     test = zeros(Height,Width);
%     counts = zeros(Height,Width);
%     
%     for x = 0:Height-1
%         for y = 1:Width-1
%             if(x == 0)
%                 test(1,y) = UF_Array(1,y+x*Height);
%                 counts(1,y) = UF_Array(2,y+x*Height);
%             else
%                 test(x,y) = UF_Array(1,y+x*Height);
%                 counts(x,y) = UF_Array(2,y+x*Height);
%             end
%         end
%     end
%     
%     test1 = diff(test,1,1);
%     test2 = diff(test,1,2);
%     test1(test1 == 480) = 0;
%     test1 = vertcat(HorzBuf,test1);
%     %test1(test1 ~= 0) = 255;
%     
%     test2(test2 == 1) = 0;
%     test2 = horzcat(VertBuf,test2);
%     %test2(test2 ~= 0) = 255;
%     
%     test12 = test1 + test2;
%     figure;
%     imshow(test12);
% end
