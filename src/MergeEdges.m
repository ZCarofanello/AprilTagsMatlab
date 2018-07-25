function Clusters = MergeEdges(Edges,Magnitude,Direction)
    %Constants to export sometime
    thetaThr = 100;
    magThr = 1200;
    
    %Get the width and height of the iamge
    width = size(Magnitude,2);
    height = size(Magnitude,1);
    
    %Reshape the Magnitude and Directions of the arrays
    tmin = ArraytoList(Direction);
    tmax = tmin;
    mmin = ArraytoList(Magnitude);
    mmax = mmin;
    
    ValidIds = [Edges(:,2) ; Edges(:,3)];
    
    %Create the unionfind vector which is pre allocated for speed
    SimpleUF = [(1:width*height)', ones(1,width*height)'];
    
    test = ismember(SimpleUF(:,1),ValidIds);
    SimpleUF(~test) = 0;
    
    for i = 1:size(Edges,1)
        ida = Edges(i,2);
        idb = Edges(i,3);
        
        ida = IgetRepresentative(SimpleUF,ida); %gets rep
        idb = IgetRepresentative(SimpleUF,idb); %gets rep
        
        if(ida == idb) %It's already connected!
            continue;
        end
        
        sza = SimpleUF(ida,2); %Get the size of tree a
        szb = SimpleUF(idb,2); %Get the size of tree b
        
        tmina = tmin(ida); tmaxa = tmax(ida); %finds the max and mins of a
        tminb = tmin(idb); tmaxb = tmax(idb); %finds the max and mins of b
        
        costa = tmaxa-tmina; %Intermediate cost value
        costb = tmaxb-tminb; %Intermediate cost value
        
        %Makes sure that the angles aren't more than +/- PI
        bshift = mod2pi((tminb+tmaxb)/2,(tmina+tmaxa)/2)-(tminb+tmaxb)/2;
        
        tminab = min(tmina, tminb+bshift); %Theta min
        tmaxab = max(tmaxa, tmaxb+bshift); %Theta max
        
        if(tmaxab-tminab > 2*pi)
            tmaxab = tminab + 2*pi;
        end
        
        mminab = min(mmin(ida), mmin(idb)); %Mag min
        mmaxab = max(mmax(ida), mmax(idb)); %Mag max
        
        costab = (tmaxab - tminab); %Intermediate cost value for a and b
        
        %Magic Values that I need to understand more :)
        Value1 = (costab <= min(costa,costb) + (thetaThr/(sza+szb)));
        Value2 = (mmaxab-mminab) <= min(mmax(ida)-mmin(ida),...
            mmax(idb)-mmin(idb)) + (magThr/(sza+szb));
        
        if(Value1 && Value2)
            [SimpleUF, idab] = IconnectNodes(SimpleUF, ida, idb,test);
            
            tmin(idab) = tminab; %Sets the minimum theta
            tmax(idab) = tmaxab; %Sets the maximum theta
            
            mmin(idab) = mminab; %Sets the minimum mag
            mmin(idab) = mmaxab; %Sets the maximum mag
        end
    end
    %Export the clusters
    Clusters = ExportClusters(SimpleUF,Magnitude, Edges);
end


% Gets the representative of the node
function root = IgetRepresentative(UFArray,NodeId)
    if(UFArray(NodeId,1) == NodeId) %If it is it's own rep return
        root = NodeId;              %No changes
    else
        root = UFArray(NodeId,1);
    end
end

% Gets the representative of the node
function [root,UpdatedArray] = getRepresentative(UFArray,NodeId)
    if(UFArray(NodeId,1) == NodeId) %If it is it's own rep return
        root = NodeId;              %No changes
    else
        root = getRepresentative(UFArray,UFArray(NodeId,1)); %Recurse
        UFArray(NodeId,1) = root; %Flatten the tree
    end
UpdatedArray = UFArray;   %Return the updated array
end

%connects and merges the two trees together
function [UFArray,root] = IconnectNodes(UFArray, aId,bId,ValidIds)

    aRoot = IgetRepresentative(UFArray,aId); %Get rep of a
    bRoot = IgetRepresentative(UFArray,bId); %Get rep of b

    if(aRoot==bRoot) %It's already connected!
        root=aRoot;  %Return the root
        return;
    end
    
    if(UFArray(aRoot,2) > UFArray(bRoot,2)) %Larger tree wins!
        %Add the sizes together
        UFArray(aRoot,2) = UFArray(aRoot,2) + UFArray(bRoot,2);
        
        UFArray(UFArray == bRoot) = aRoot;
        
        root=aRoot; %Return the new root
        return;
    else
        %Add the sizes together
        UFArray(bRoot,2) = UFArray(aRoot,2) + UFArray(bRoot,2);
        
        UFArray(UFArray == aRoot) = bRoot;
        
        root=bRoot; %Return the new root
        return;
    end
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
    
    if(UFArray(aRoot,2) > UFArray(bRoot,2)) %Larger tree wins!
        NewUFArray(bRoot,1) = aRoot; %Set the new root
        
        %Add the sizes together
        NewUFArray(aRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2);
        
        root=aRoot; %Return the new root
        return;
     else
        NewUFArray(aRoot,1) = bRoot; %Set the new root
        
        %Add the sizes together
        NewUFArray(bRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2);
        
        root=bRoot; %Return the new root
        return;
    end
end

function FixedTree = FlattenTree(OldRoot,NewRoot,UF_Array)
Children_Idx = FindChildren(UF_Array,OldRoot); %Gets all children from root
UF_Array(Children_Idx,1) = NewRoot;            %Shortcuts all children
FixedTree = UF_Array;
end

function Nodes = FindChildren(UF_Array, ParentIdx)
Nodes = find(UF_Array(:,1) == ParentIdx); %Finds all the children Idx
end

function longArray = ArraytoList(Array)
Width = size(Array,2);
Height  = size(Array,1);

longArray = zeros(1,Width*Height);
for i = 1:Height
    StartIdx = ((i-1) * Width)+1;
    EndIdx   = (StartIdx + Width)-1;
    longArray(1,StartIdx:EndIdx) = Array(i,:);
end
end

%Formatting the clusters as a list with points to make it easier later
function ClusterList = ExportClusters(UF_Array,Magnitude,Edges)
    %Need to export these constants
    MinCluster = 4;

    %find clusters that have more than the MinSeg
    Valid_Clusters = UF_Array((UF_Array(:,2) >= MinCluster),1);
    
    %Extra check to flatten tree (not necessary)
%     for k = 1:length(Valid_Clusters)
%         root = getRepresentative(UF_Array,Valid_Clusters(k));
%         UF_Array = FlattenTree(Valid_Clusters(k),root,UF_Array);
%     end

    %Create a logical array for faster indexing / display
    logical_arr = ismember(UF_Array(:,1),Valid_Clusters);
    
    ClusterList = [];  %Empty matrix for clusters
    FirstEntry = true; %Bool to make sure we don't miss the first entry
    
    for i = 1:size(Edges,1)-1 %loops through all the edges
        if(logical_arr(Edges(i,2))) %Is the edge a part of valid cluster
            
            EdgeCluster = UF_Array(Edges(i,3),1); %Gets cluster #
            EdgeMag = Magnitude(Edges(i,3));      %Gets magnitude
            EdgeX = Edges(i,4);                   %Gets X coord
            EdgeY = Edges(i,5);                   %Gets Y coord
            if(FirstEntry)
                ClusterList = [EdgeX,EdgeY,EdgeMag,EdgeCluster];
                FirstEntry = false;
            else
                ClusterList=[ClusterList;EdgeX,EdgeY,EdgeMag,EdgeCluster];
            end
        end
    end

    ClusterList = sortrows(ClusterList,4);
end